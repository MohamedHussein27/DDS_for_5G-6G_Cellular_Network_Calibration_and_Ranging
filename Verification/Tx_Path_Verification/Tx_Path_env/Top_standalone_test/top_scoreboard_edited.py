"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
 
    Module: tx_scoreboard.py
 
    Description:
        The TX scoreboard receives paired transactions from the monitor.
        The monitor captures every clock cycle and packs the DUT inputs
        and outputs into a top_item transaction.
 
        Unlike the IFFT scoreboard which accumulates N samples before
        computing, here the golden model only needs (FTW_start, FTW_step,
        cycles) once to produce the full output vectors. The scoreboard
        calls the golden model as soon as dds_enable first asserts,
        then uses a simple index to compare each DUT output sample
        against the pre-computed reference vector as tx_valid pulses.
"""
 
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import logging
import numpy as np
 
from top_seq_item import *
 
# Import your TX golden model
from tx_golden_model import *
 
# ─────────────────────────────────────────────────────────────────────────────
# RTL parameters  (mirror Tx_path_fixed.py exactly)
# ─────────────────────────────────────────────────────────────────────────────
WL        = 16      # word length (signed fixed-point)
N         = 4096    # FFT / IFFT / MUX frame size
FL        = 8       # fractional bits  (Q8.8 at TX output)
TOLERANCE = 4       # ±LSBs allowed between DUT and golden model
 
 
# ─────────────────────────────────────────────────────────────────────────────
# Fixed-point helpers  (identical to IFFT scoreboard helpers)
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
    raw = max(-(1 << (WL - 1)), min((1 << (WL - 1)) - 1, raw))
    return raw
 
 
# ─────────────────────────────────────────────────────────────────────────────
# Scoreboard
# ─────────────────────────────────────────────────────────────────────────────
class top_scoreboard(uvm_scoreboard):
    """
    Scoreboard for the full TX path (TX_TOP).
 
    Inputs captured by the monitor in top_item:
        rst_n, dds_ready_flag, FTW_start, FTW_step, cycles
 
    Outputs captured by the monitor in top_item:
        tx_valid, tx_out_real, tx_out_imag
    """
 
    def __init__(self, name, parent):
        super().__init__(name, parent)
 
    # ── build_phase ───────────────────────────────────────────────────────────
    def build_phase(self):
        # top analysis export and FIFO for receiving transactions from the monitor
        self.sb_export = uvm_analysis_export("sb_export", self)
        self.sb_fifo   = uvm_tlm_analysis_fifo("sb_fifo", self)
        self.sb_export = self.sb_fifo.analysis_export

        # DDS analysis port (to receive FTW and cycle info for golden model)
        self.dds_export = uvm_analysis_port("dds_export", self)
        self.dds_fifo   = uvm_tlm_analysis_fifo("dds_fifo", self)
        self.dds_export = self.dds_fifo.analysis_export
 
        # Counters
        self.correct_real = 0
        self.error_real   = 0
        self.correct_imag = 0
        self.error_imag   = 0
        self.reset_cycles = 0
 
        # Golden reference vectors (filled once on first dds_enable)
        self._ref_real  = []    # list of N raw signed integers
        self._ref_imag  = []    # list of N raw signed integers
        self._golden_ready    = False   # True after golden model has been called
        self._out_idx         = 0       # which sample we are comparing next
 
    # ── run_phase ─────────────────────────────────────────────────────────────
    async def run_phase(self):
        while True:
            item = await self.sb_fifo.get()
            self._process(item)
 
    # ── _process ──────────────────────────────────────────────────────────────
    def _process(self, item):
 
        global i
 
        # ── Reset: flush state ────────────────────────────────────────────────
        if not item.rst_n:
            self.reset_cycles += 1
            self._ref_real     = []
            self._ref_imag     = []
            self._golden_ready = False
            self._out_idx      = 0
            self.logger.info("Reset detected, golden model flushed.")
            i = 0
            return
 
        # ── Call golden model once on the first dds_enable pulse ──────────────
        # The golden model only needs (FTW_start, FTW_step, cycles) to produce
        # the full N-sample output vectors – no accumulation needed.
        if item.dds_ready_flag and not self._golden_ready:
            self.logger.info(
                f"Calling golden model: "
                f"FTW_start={item.FTW_start}, FTW_step={item.FTW_step}"
            )
 
            # tx_golden_model returns two vectors: real and imag (raw signed ints)
            ref_real_vec, ref_imag_vec = run_tx_top_pipeline(
                item.FTW_start, item.FTW_step, N_cycles=4096,Fs=491.52e6,
                ofdm_re_file="ofdm_data_re.hex",
                ofdm_im_file="ofdm_data_im.hex",
            )
 
            self._ref_real     = list(ref_real_vec)
            self._ref_imag     = list(ref_imag_vec)
            self._golden_ready = True
            self._out_idx      = 0
 
            # Debug: print first and last reference samples
            print(f"Golden model output sample #0:    re={self._ref_real[0]}    im={self._ref_imag[0]}")
            print(f"Golden model output sample #4095: re={self._ref_real[4095]} im={self._ref_imag[4095]}")
 
        # ── Debug heartbeat / boundary prints ────────────────────────────────
        if (i > 4090 and i < 4100) or (i > 12275 and i < 12285):
            self.logger.info(
                f"TX scoreboard sample #{i}: "
                f"| DUT output: re={_sign16(item.tx_out_real)}, "
                f"im={_sign16(item.tx_out_imag)}, "
                f"dds_ready_flag={item.dds_ready_flag}, tx_valid={item.tx_valid}"
            )
 
        if i % 1000 == 0:
            self.logger.info(f"--- SIMULATION HEARTBEAT: Processing cycle {i} ---")
 
        # ── Compare DUT output with reference vector ───────────────────────────
        if item.tx_valid:
            if not self._golden_ready:
                self.logger.warning(
                    "DUT asserted tx_valid but golden model has not been called yet. "
                    "Possible timing issue – skipping comparison."
                )
            elif self._out_idx >= N:
                self.logger.warning(
                    f"DUT asserted tx_valid beyond expected frame size N={N} "
                    f"(out_idx={self._out_idx}) – skipping comparison."
                )
            else:
                dut_real = _sign16(item.tx_out_real)
                dut_imag = _sign16(item.tx_out_imag)
                ref_real = self._ref_real[self._out_idx]
                ref_imag = self._ref_imag[self._out_idx]
 
                # Debug print at frame boundaries
                if (i > 4090 and i < 4100):
                    self.logger.info(
                        f"Comparing sample #{self._out_idx}: "
                        f"DUT(re={dut_real}, im={dut_imag}) vs "
                        f"REF(re={ref_real}, im={ref_imag})"
                    )
 
                self._compare_signal("tx_out_real", dut_real, ref_real)
                self._compare_signal("tx_out_imag", dut_imag, ref_imag)
                self._out_idx += 1
 
        # Update global cycle counter
        i += 1
 
    # ── _compare_signal ───────────────────────────────────────────────────────
    def _compare_signal(self, signal_name, dut_val, ref_val):
        """Compare one signal with ±TOLERANCE LSB window."""
        diff = abs(dut_val - ref_val)
        if diff <= TOLERANCE:
            if signal_name == "tx_out_real":
                self.correct_real += 1
            else:
                self.correct_imag += 1
        else:
            if signal_name == "tx_out_real":
                self.error_real += 1
            else:
                self.error_imag += 1
            self.logger.error(
                f"MISMATCH {signal_name} "
                f"(error #{self.error_real if signal_name == 'tx_out_real' else self.error_imag}) "
                f"| DUT={dut_val}  REF={ref_val}  diff={diff}  tolerance={TOLERANCE}"
            )
 
    # ── report_phase ──────────────────────────────────────────────────────────
    def report_phase(self):
        total_real = self.correct_real + self.error_real
        total_imag = self.correct_imag + self.error_imag
        pass_real  = (self.error_real == 0)
        pass_imag  = (self.error_imag == 0)
        overall    = "PASS ✓" if (pass_real and pass_imag) else "FAIL ✗"
 
        self.logger.info("╔══════════════════════════════════════════════════════╗")
        self.logger.info("║              TX SCOREBOARD REPORT                    ║")
        self.logger.info("╠══════════════════════════════════════════════════════╣")
        self.logger.info(f"║  Reset cycles detected        : {self.reset_cycles:<22} ║")
        self.logger.info("║ ───────────────────────────────────────────────────  ║")
        self.logger.info(f"║  tx_out_real correct / total  : {self.correct_real}/{total_real:<20} ║")
        self.logger.info(f"║  tx_out_real errors           : {self.error_real:<22} ║")
        self.logger.info("║ ───────────────────────────────────────────────────  ║")
        self.logger.info(f"║  tx_out_imag correct / total  : {self.correct_imag}/{total_imag:<20} ║")
        self.logger.info(f"║  tx_out_imag errors           : {self.error_imag:<22} ║")
        self.logger.info("║ ───────────────────────────────────────────────────  ║")
        self.logger.info(f"║  OVERALL RESULT               : {overall:<22} ║")
        self.logger.info("╚══════════════════════════════════════════════════════╝")
 
        # Hard-fail the test if there were any mismatches
        if not (pass_real and pass_imag):
            self.logger.critical(
                f"Scoreboard detected {self.error_real} tx_out_real error(s) "
                f"and {self.error_imag} tx_out_imag error(s)."
            )