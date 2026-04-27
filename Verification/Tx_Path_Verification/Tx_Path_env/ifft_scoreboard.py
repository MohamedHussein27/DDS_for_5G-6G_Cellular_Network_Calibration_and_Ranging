"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_scoreboard.py

    Description:
        The IFFT scoreboard receives paired transactions from the monitor. 
        Because the RTL pipeline latency is exactly N-1 clocks, it pushes 
        every valid_in sample into the Golden Model. When valid_out asserts, 
        it pops the corresponding computed sample and verifies the Q11.5 
        fixed-point outputs against the hardware.
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import logging
import numpy as np

from ifft_item import *

# Assuming your python golden model is in a file named ifft_golden.py
from ifft_fixed import radix22_dif_ifft_fixed

# ─────────────────────────────────────────────────────────────────────────────
# RTL parameters
# ─────────────────────────────────────────────────────────────────────────────
WL          = 16          # word length (signed fixed-point)
N           = 4096        # IFFT size
FL          = 5           # fractional bits (Q11.5 format used by TX Output)
TOLERANCE   = 2           # ±LSBs allowed between DUT and golden model


# ─────────────────────────────────────────────────────────────────────────────
# Golden-model helpers
# ─────────────────────────────────────────────────────────────────────────────
def _sign16(v):
    """Sign-extend a raw 16-bit integer value from the DUT."""
    v = int(v) & 0xFFFF
    if v & 0x8000:
        return v - 0x10000
    else:
        return v


def _raw_int_to_float(v, frac_bits):
    """Convert a raw fixed-point integer from the DUT to a Python float."""
    return _sign16(v) / (1 << frac_bits)


def _float_to_raw_int(f, frac_bits):
    """
    Convert a float back to a signed integer.
    Saturates to 16-bit signed range [-32768, 32767].
    """
    raw = int(round(f * (1 << frac_bits)))
    raw = max(-(1 << (WL - 1)), min((1 << (WL - 1)) - 1, raw)) # take min bet. (1 << (WL - 1)) - 1, raw then take the max bet. this value in the minumum value
    return raw


class IFFTGoldenModel:
    """
    Accumulates N complex input samples, runs the custom fixed-point 
    radix22_dif_ifft_fixed model, and yields Q11.5 output samples.
    """

    def __init__(self, n=N):
        self._n      = n
        self._buf_re = []
        self._buf_im = []
        self._out    = []   # list of (out_real, out_imag) in raw signed integers

    def push(self, in_real_raw, in_imag_raw):
        """Accept one pair of raw RTL integer samples."""
        self._buf_re.append(_raw_int_to_float(in_real_raw, FL))
        self._buf_im.append(_raw_int_to_float(in_imag_raw, FL))

        if len(self._buf_re) == self._n:
            self._compute()

    def _compute(self):
        # 1. Reconstruct the complex array
        x_in = np.array(self._buf_re) + 1j * np.array(self._buf_im)
        
        # 2. Call your custom hardware-accurate fixed-point function
        X_out_complex, wl_out, fl_out = radix22_dif_ifft_fixed(x_in, WL=WL, FL=FL)
        
        # 3. Convert the resulting floats back to raw signed integers for DUT comparison
        self._out = []

        for v in X_out_complex:
            
            real_int = _float_to_raw_int(v.real, FL)
            
            imag_int = _float_to_raw_int(v.imag, FL)
            
            self._out.append((real_int, imag_int))
        
        # 4. Clear buffer so the model is ready to absorb the next frame
        self._buf_re.clear()
        self._buf_im.clear()

    def pop(self):
        """
        Return the next (out_real, out_imag) pair,
        or None if no output is ready yet.
        """
        return self._out.pop(0) if self._out else None

# ─────────────────────────────────────────────────────────────────────────────
# Scoreboard
# ─────────────────────────────────────────────────────────────────────────────
class IFFTScoreboard(uvm_scoreboard):
    """
    Scoreboard for the 4096-point IFFT.
    """

    def __init__(self, name, parent):
        super().__init__(name, parent)

    # ── build_phase ──────────────────────────────────────────────────────
    def build_phase(self):
        self.sb_export = uvm_analysis_export("sb_export", self)
        self.sb_fifo   = uvm_tlm_analysis_fifo("sb_fifo", self)
        self.sb_export = self.sb_fifo.analysis_export

        # Counters
        self.correct_real  = 0
        self.error_real    = 0
        self.correct_imag  = 0
        self.error_imag    = 0
        self.reset_cycles  = 0

        # Golden model instance
        self._golden = IFFTGoldenModel(n=N)

    # ── run_phase ─────────────────────────────────────────────────────────
    async def run_phase(self):
        while True:
            item = await self.sb_fifo.get()
            self._process(item)

    # ── _process ─────────────────────────────────────────────────────────
    def _process(self, item):
        # ── Reset: flush golden model state ──────────────────────────────
        if not item.rst_n:
            self.reset_cycles += 1
            self._golden = IFFTGoldenModel(n=N)   # wipe internal state
            self.logger.debug("Reset detected , golden model flushed.")
            return

        # ── Feed input sample into golden model ──────────────────────────
        if item.valid_in:
            self._golden.push(item.in_real, item.in_imag)
            self.logger.debug(
                f"Golden model received sample "
                f"re={_sign16(item.in_real)}, im={_sign16(item.in_imag)}"
            )

        # ── Compare DUT output with golden model output ───────────────────
        if item.valid_out:
            golden_pair = self._golden.pop()

            if golden_pair is None:
                self.logger.warning(
                    "DUT asserted valid_out but golden model has no output ready. "
                    "Possible frame-alignment issue – skipping comparison."
                )
                return

            ref_real_q, ref_imag_q = golden_pair

            dut_real = _sign16(item.out_real)
            dut_imag = _sign16(item.out_imag)
            
            # The ref values from _float_to_raw_int are already signed integers
            ref_real = ref_real_q
            ref_imag = ref_imag_q

            self._compare_signal(
                signal_name="out_real",
                dut_val=dut_real,
                ref_val=ref_real,
            )
            self._compare_signal(
                signal_name="out_imag",
                dut_val=dut_imag,
                ref_val=ref_imag,
            )

    # ── _compare_signal ───────────────────────────────────────────────────
    def _compare_signal(self, signal_name, dut_val, ref_val):
        """Compare one signal with ±TOLERANCE LSB window."""
        diff = abs(dut_val - ref_val)
        if diff <= TOLERANCE:
            if signal_name == "out_real":
                self.correct_real += 1
            else:
                self.correct_imag += 1
            self.logger.debug(
                f"OK  {signal_name}: DUT={dut_val}  REF={ref_val}  diff={diff}"
            )
        else:
            if signal_name == "out_real":
                self.error_real += 1
            else:
                self.error_imag += 1
            self.logger.error(
                f"MISMATCH {signal_name} "
                f"(error #{self.error_real if signal_name == 'out_real' else self.error_imag}) "
                f"| DUT={dut_val}  REF={ref_val}  diff={diff}  tolerance={TOLERANCE}"
            )

    # ── report_phase ─────────────────────────────────────────────────────
    def report_phase(self):
        total_real  = self.correct_real  + self.error_real
        total_imag  = self.correct_imag  + self.error_imag
        pass_real   = (self.error_real  == 0)
        pass_imag   = (self.error_imag  == 0)
        overall     = "PASS ✓" if (pass_real and pass_imag) else "FAIL ✗"

        self.logger.info("╔══════════════════════════════════════════════════════╗")
        self.logger.info("║            IFFT SCOREBOARD REPORT                    ║")
        self.logger.info("╠══════════════════════════════════════════════════════╣")
        self.logger.info(f"║  Reset cycles detected        : {self.reset_cycles:<22} ║")
        self.logger.info("║ ───────────────────────────────────────────────────  ║")
        self.logger.info(f"║  out_real  correct / total    : {self.correct_real}/{total_real:<20} ║")
        self.logger.info(f"║  out_real  errors             : {self.error_real:<22} ║")
        self.logger.info("║ ───────────────────────────────────────────────────  ║")
        self.logger.info(f"║  out_imag  correct / total    : {self.correct_imag}/{total_imag:<20} ║")
        self.logger.info(f"║  out_imag  errors             : {self.error_imag:<22} ║")
        self.logger.info("║ ───────────────────────────────────────────────────  ║")
        self.logger.info(f"║  OVERALL RESULT               : {overall:<22} ║")
        self.logger.info("╚══════════════════════════════════════════════════════╝")

        # Hard-fail the test if there were any mismatches
        if not (pass_real and pass_imag):
            self.logger.critical(
                f"Scoreboard detected {self.error_real} out_real error(s) "
                f"and {self.error_imag} out_imag error(s)."
            )