"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: fft_scoreboard.py

    Description:
        The FFT scoreboard receives paired transactions from the monitor. 
        It pushes valid_in samples into the bit-true Golden Model. 
        When valid_out asserts, it pops the computed sample and verifies 
        the raw fixed-point integer outputs against the hardware.
"""
import cocotb
from cocotb.triggers import *
from pyuvm import *
import logging
import numpy as np
import os

from fft_seq_item import fft_item
from fft_golden_model import radix2_dif_fft_fixed

# ─────────────────────────────────────────────────────────────────────────────
# RTL parameters
# ─────────────────────────────────────────────────────────────────────────────
WL          = 16          # word length (signed fixed-point)
N           = 4096        # FFT size
TOLERANCE   = 4           # ±LSBs allowed between DUT and golden model

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

class FFTGoldenModel:
    """
    Accumulates N complex input samples, runs the custom fixed-point 
    radix2_dif_fft_fixed model, and queues raw integer output samples.
    """
    def __init__(self, n=N):
        self._n      = n
        self._buf_re = []
        self._buf_im = []
        self._out    = []   # list of (out_real, out_imag) in raw signed integers
        self.frames  = 0    # Track how many frames we've computed

    def push(self, in_real_raw, in_imag_raw):
        """Accept one pair of raw RTL integer samples."""
        # 1-Cycle Register Delay handling (if RTL still requires it!)
        # If your RTL is perfectly aligned to valid_in now, you can remove this comment block.
        self._buf_re.append(_sign16(in_real_raw))
        self._buf_im.append(_sign16(in_imag_raw))

        if len(self._buf_re) == self._n:
            self.frames += 1
            logging.info(f"Golden model buffer full. Computing FFT Frame {self.frames}...")
            self._compute()

    def _compute(self):
        # 1. Reconstruct the complex integer array
        x_in = np.array(self._buf_re) + 1j * np.array(self._buf_im)
        
        # 2. Call the hardware-accurate fixed-point function (Expects raw integers)
        X_out_complex = radix2_dif_fft_fixed(x_in, WL=WL, is_ifft=False)

        # 3. Store raw integers directly to output queue
        self._out = []
        for v in X_out_complex:
            self._out.append((int(v.real), int(v.imag)))
        
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
class FFTScoreboard(uvm_scoreboard):
    """
    Scoreboard for the 4096-point FFT.
    """
    def __init__(self, name, parent):
        super().__init__(name, parent)
    def build_phase(self):
        self.sb_fifo   = uvm_tlm_analysis_fifo("sb_fifo", self)
        self.sb_export = self.sb_fifo.analysis_export

        # Counters
        self.correct_real  = 0
        self.error_real    = 0
        self.correct_imag  = 0
        self.error_imag    = 0
        self.reset_cycles  = 0
        self.output_idx    = 0

        # Golden model instance
        self._golden = FFTGoldenModel(n=N)

        # ─────────────────────────────────────────────────────────
        # EXTERNAL FILE DUMP SETUP
        # ─────────────────────────────────────────────────────────
        self.f_rtl = open("rtl_out.csv", "w")
        self.f_ref = open("python_out.csv", "w")
        
        # Write the CSV headers
        self.f_rtl.write("RTL_Real,RTL_Imag\n")
        self.f_ref.write("REF_Real,REF_Imag\n")

    async def run_phase(self):
        while True:
            item = await self.sb_fifo.get()
            self._process(item)

    def _process(self, item):
        # ── Reset: flush golden model state ──────────────────────────────
        if not item.rst_n:
            self.reset_cycles += 1
            self._golden = FFTGoldenModel(n=N)   # wipe internal state
            self.output_idx = 0
            self.logger.info("Reset detected, golden model flushed.")
            return

        # ── Feed input sample into golden model ──────────────────────────
        if item.valid_in:
            self._golden.push(item.in_real, item.in_imag)

        # ── Compare DUT output with golden model output ───────────────────
        if item.valid_out:
            golden_pair = self._golden.pop()

            if golden_pair is None:
                self.logger.warning(
                    "DUT asserted valid_out but golden model has no output ready. "
                    "Possible frame-alignment issue – skipping comparison."
                )
                return

            ref_real, ref_imag = golden_pair
            dut_real = _sign16(item.out_real)
            dut_imag = _sign16(item.out_imag)

            # Debug print for first 5 samples of every frame
            # ─────────────────────────────────────────────────────────
            # DUMP ALL 4096 SAMPLES TO TERMINAL
            # ─────────────────────────────────────────────────────────
            self.logger.info(
                f"Sample #{self.output_idx:04d} | "
                f"REF: re={ref_real:6d}, im={ref_imag:6d} | "
                f"DUT: re={dut_real:6d}, im={dut_imag:6d}"
            )
            # ─────────────────────────────────────────────────────────
            # WRITE TO EXTERNAL FILES
            # ─────────────────────────────────────────────────────────
            self.f_rtl.write(f"{dut_real},{dut_imag}\n")
            self.f_ref.write(f"{ref_real},{ref_imag}\n")

            self._compare_signal("out_real", dut_real, ref_real)
            self._compare_signal("out_imag", dut_imag, ref_imag)

            self.output_idx += 1
            if self.output_idx == N:
                self.output_idx = 0

    def _compare_signal(self, signal_name, dut_val, ref_val):
        """Compare one signal with ±TOLERANCE LSB window."""
        diff = abs(dut_val - ref_val)
        if diff <= TOLERANCE:
            if signal_name == "out_real":
                self.correct_real += 1
            else:
                self.correct_imag += 1
        else:
            if signal_name == "out_real":
                self.error_real += 1
            else:
                self.error_imag += 1
            self.logger.error(
                f"MISMATCH {signal_name} "
                f"(error #{self.error_real if signal_name == 'out_real' else self.error_imag}) "
                f"| Index={self.output_idx} | DUT={dut_val}  REF={ref_val}  diff={diff}"
            )

    def report_phase(self):
        total_real  = self.correct_real  + self.error_real
        total_imag  = self.correct_imag  + self.error_imag
        pass_real   = (self.error_real  == 0)
        pass_imag   = (self.error_imag  == 0)
        overall     = "PASS ✓" if (pass_real and pass_imag) else "FAIL ✗"

        self.logger.info("")
        self.logger.info("=== FFT SCOREBOARD REPORT ===")
        self.logger.info(f" Reset Cycles : {self.reset_cycles}")
        self.logger.info(f" Real Output  : {self.correct_real}/{total_real} Correct | {self.error_real} Errors")
        self.logger.info(f" Imag Output  : {self.correct_imag}/{total_imag} Correct | {self.error_imag} Errors")
        self.logger.info(f" OVERALL      : {overall}")
        self.logger.info("=============================")
        self.logger.info("")

        if not (pass_real and pass_imag):
            self.logger.critical(
                f"Scoreboard detected {self.error_real} out_real error(s) "
                f"and {self.error_imag} out_imag error(s)."
            )