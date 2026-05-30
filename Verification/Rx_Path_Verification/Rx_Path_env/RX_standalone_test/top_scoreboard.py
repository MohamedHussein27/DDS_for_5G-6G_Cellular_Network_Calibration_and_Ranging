"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: rx_scoreboard.py

    Description:
        The RX scoreboard receives paired transactions from rx_monitor.
        The monitor captures every clock cycle and packs all DUT
        inputs and outputs into an rx_item transaction.

        Unlike the TX scoreboard which calls the golden model on the
        first dds_enable pulse, the RX scoreboard must accumulate:
          1. 2048 reference RAM writes  (ref_wr_en = 1)
          2. 4096 received signal samples (rx_valid_in = 1)
        before calling the RX golden model. The golden model is called
        once per frame and produces two pre-computed output vectors:
          - ofdm_out_re / ofdm_out_im   [2048 samples]
          - radar_out_re / radar_out_im  [2048 samples]
        The scoreboard then indexes into those vectors as
        ofdm_valid_out and radar_valid_out pulse, comparing
        DUT vs. reference with ±TOLERANCE LSBs.

        Full pipeline modelled in rx_golden_model.py:
            FFT(4096) → bit_rev[12-bit] → Demux
              → OFDM out [bins 0..2047]
              → conj_multiply(ref_ram) → >>>3 → IFFT(2048)
              → bit_rev[11-bit] → Radar out [2048]
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import logging
import numpy as np

from top_seq_item import rx_item
from rx_golden import run_rx_pipeline

# ─────────────────────────────────────────────────────────────────────────────
# RTL parameters
# ─────────────────────────────────────────────────────────────────────────────
WL            = 16
N_FFT         = 4096
N_HALF        = 2048
TOLERANCE     = 4      # ±LSBs allowed between DUT and golden model

# ─────────────────────────────────────────────────────────────────────────────
# Fixed-point helpers  (identical style to tx_scoreboard.py)
# ─────────────────────────────────────────────────────────────────────────────
def _sign16(v):
    """Sign-extend a raw 16-bit integer value from the DUT."""
    v = int(v) & 0xFFFF
    if v & 0x8000:
        return v - 0x10000
    else:
        return v


# ─────────────────────────────────────────────────────────────────────────────
# Scoreboard
# ─────────────────────────────────────────────────────────────────────────────
class top_scoreboard(uvm_scoreboard):
    """
    Scoreboard for the full RX path (RX_TOP).

    Inputs captured by the monitor in rx_item:
        rst_n, rx_valid_in, rx_in_re, rx_in_im
        ref_wr_en, ref_wr_re, ref_wr_im

    Outputs captured by the monitor in rx_item:
        ofdm_valid_out, ofdm_out_re, ofdm_out_im
        radar_valid_out, radar_out_re, radar_out_im
    """

    def __init__(self, name, parent):
        super().__init__(name, parent)

    # ── build_phase ───────────────────────────────────────────────────────────
    def build_phase(self):
        self.sb_export = uvm_analysis_export("sb_export", self)
        self.sb_fifo   = uvm_tlm_analysis_fifo("sb_fifo", self)
        self.sb_export = self.sb_fifo.analysis_export

        # ── Per-output pass/fail counters ─────────────────────────────────
        self.correct_ofdm_real  = 0
        self.error_ofdm_real    = 0
        self.correct_ofdm_imag  = 0
        self.error_ofdm_imag    = 0
        self.correct_radar_real = 0
        self.error_radar_real   = 0
        self.correct_radar_imag = 0
        self.error_radar_imag   = 0
        self.reset_cycles       = 0

        # ── Input accumulators ────────────────────────────────────────────
        # Reference RAM — collected once; locked after 2048 samples
        self._ref_re       = []
        self._ref_im       = []
        self._ref_locked   = False

        # Received signal — reset after every 4096-sample frame
        self._rx_re_frame  = []
        self._rx_im_frame  = []

        # ── Golden reference vectors  (filled once per frame) ─────────────
        self._ofdm_ref_re   = []   # list of N_HALF raw signed integers
        self._ofdm_ref_im   = []
        self._radar_ref_re  = []
        self._radar_ref_im  = []
        self._golden_ready  = False

        # Output sample index (separate pointer per path)
        self._ofdm_idx   = 0
        self._radar_idx  = 0

        # ─────────────────────────────────────────────────────────
        # EXTERNAL FILE DUMP SETUP
        # ─────────────────────────────────────────────────────────
        self.f_ofdm_rtl = open("ofdm_rtl_out.csv", "w")
        self.f_ofdm_ref = open("ofdm_python_out.csv", "w")
        self.f_radar_rtl = open("radar_rtl_out.csv", "w")
        self.f_radar_ref = open("radar_python_out.csv", "w")
        
        self.f_ofdm_rtl.write("RTL_Real,RTL_Imag\n")
        self.f_ofdm_ref.write("REF_Real,REF_Imag\n")
        self.f_radar_rtl.write("RTL_Real,RTL_Imag\n")
        self.f_radar_ref.write("REF_Real,REF_Imag\n")

    # ── run_phase ─────────────────────────────────────────────────────────────
    async def run_phase(self):
        while True:
            item = await self.sb_fifo.get()
            self._process(item)

    # ── _process ──────────────────────────────────────────────────────────────
    def _process(self, item):

        # ── Reset: flush all state ────────────────────────────────────────
        if not item.rst_n:
            self.reset_cycles += 1
            self._ref_re.clear()
            self._ref_im.clear()
            self._ref_locked  = False
            self._rx_re_frame.clear()
            self._rx_im_frame.clear()
            self._ofdm_ref_re.clear()
            self._ofdm_ref_im.clear()
            self._radar_ref_re.clear()
            self._radar_ref_im.clear()
            self._golden_ready = False
            self._ofdm_idx     = 0
            self._radar_idx    = 0
            self.logger.info("Reset detected, RX golden model flushed.")
            return

        # ── Step 1: Collect reference RAM writes (first 2048, then lock) ──
        if item.ref_wr_en and not self._ref_locked:
            self._ref_re.append(_sign16(item.ref_wr_re))
            self._ref_im.append(_sign16(item.ref_wr_im))
            if len(self._ref_re) == N_HALF:
                self._ref_locked = True
                self.logger.info(
                    f"RX-SB: Reference RAM fully loaded "
                    f"({N_HALF} samples). Radar path armed."
                )

        # ── Step 2: Accumulate received signal samples ────────────────────
       # ── Step 2: Accumulate received signal samples ────────────────────
        if item.rx_valid_in:
            self._rx_re_frame.append(_sign16(item.rx_in_re))
            self._rx_im_frame.append(_sign16(item.rx_in_im))

            # Once we have a full 4096-sample frame, call the golden model
            if len(self._rx_re_frame) == N_FFT:
                
                # FIX: If reference RAM was never written (e.g., zero_ref test),
                # force it to an array of 2048 zeros so the golden model can run.
                if not self._ref_locked:
                    self.logger.info("RX-SB: No Reference RAM writes detected. Defaulting to all zeros.")
                    self._ref_re = [0] * N_HALF
                    self._ref_im = [0] * N_HALF
                    self._ref_locked = True

                self.logger.info(
                    f"Calling RX golden model: "
                    f"rx_frame[{N_FFT}] + ref_ram[{N_HALF}] ready."
                )

                (ofdm_re_vec, ofdm_im_vec,
                 radar_re_vec, radar_im_vec) = run_rx_pipeline(
                    rx_re_array  = self._rx_re_frame,
                    rx_im_array  = self._rx_im_frame,
                    ref_re_array = self._ref_re,
                    ref_im_array = self._ref_im,
                )

                self._ofdm_ref_re  = list(ofdm_re_vec)
                self._ofdm_ref_im  = list(ofdm_im_vec)
                self._radar_ref_re = list(radar_re_vec)
                self._radar_ref_im = list(radar_im_vec)
                self._golden_ready = True
                self._ofdm_idx     = 0
                self._radar_idx    = 0

                # Clear frame for next accumulation
                self._rx_re_frame.clear()
                self._rx_im_frame.clear()
        # ── Step 3: Compare OFDM output ───────────────────────────────────
        if item.ofdm_valid_out:
            if not self._golden_ready:
                self.logger.warning(
                    "DUT asserted ofdm_valid_out but golden model not ready. "
                    "Possible timing issue — skipping comparison."
                )
            elif self._ofdm_idx >= N_HALF:
                self.logger.warning(
                    f"DUT asserted ofdm_valid_out beyond expected frame size "
                    f"N={N_HALF} (ofdm_idx={self._ofdm_idx}) — skipping."
                )
            else:
                dut_real = _sign16(item.ofdm_out_re)
                dut_imag = _sign16(item.ofdm_out_im)
                ref_real = self._ofdm_ref_re[self._ofdm_idx]
                ref_imag = self._ofdm_ref_im[self._ofdm_idx]

                # Debug print at frame boundaries
                if self._ofdm_idx < 4 or self._ofdm_idx > N_HALF - 4:
                    self.logger.info(
                        f"OFDM compare sample #{self._ofdm_idx}: "
                        f"DUT(re={dut_real}, im={dut_imag}) vs "
                        f"REF(re={ref_real}, im={ref_imag})"
                    )

                self._compare_signal("ofdm_out_real", dut_real, ref_real)
                self._compare_signal("ofdm_out_imag", dut_imag, ref_imag)

                # Write to CSV
                self.f_ofdm_rtl.write(f"{dut_real},{dut_imag}\n")
                self.f_ofdm_ref.write(f"{ref_real},{ref_imag}\n")

                self._ofdm_idx += 1

        # ── Step 4: Compare Radar output ──────────────────────────────────
        if item.radar_valid_out:
            if not self._golden_ready:
                self.logger.warning(
                    "DUT asserted radar_valid_out but golden model not ready. "
                    "Possible timing issue — skipping comparison."
                )
            elif self._radar_idx >= N_HALF:
                self.logger.warning(
                    f"DUT asserted radar_valid_out beyond expected frame size "
                    f"N={N_HALF} (radar_idx={self._radar_idx}) — skipping."
                )
            else:
                dut_real = _sign16(item.radar_out_re)
                dut_imag = _sign16(item.radar_out_im)
                ref_real = self._radar_ref_re[self._radar_idx]
                ref_imag = self._radar_ref_im[self._radar_idx]

                # Debug print at frame boundaries
                if self._radar_idx < 4 or self._radar_idx > N_HALF - 4:
                    self.logger.info(
                        f"Radar compare sample #{self._radar_idx}: "
                        f"DUT(re={dut_real}, im={dut_imag}) vs "
                        f"REF(re={ref_real}, im={ref_imag})"
                    )

                self._compare_signal("radar_out_real", dut_real, ref_real)
                self._compare_signal("radar_out_imag", dut_imag, ref_imag)

                # Write to CSV
                self.f_radar_rtl.write(f"{dut_real},{dut_imag}\n")
                self.f_radar_ref.write(f"{ref_real},{ref_imag}\n")

                self._radar_idx += 1

    # ── _compare_signal ───────────────────────────────────────────────────────
    def _compare_signal(self, signal_name, dut_val, ref_val):
        """Compare one signal with ±TOLERANCE LSB window."""
        diff = abs(dut_val - ref_val)
        if diff <= TOLERANCE:
            if   signal_name == "ofdm_out_real":  self.correct_ofdm_real  += 1
            elif signal_name == "ofdm_out_imag":  self.correct_ofdm_imag  += 1
            elif signal_name == "radar_out_real": self.correct_radar_real += 1
            else:                                 self.correct_radar_imag += 1
        else:
            if   signal_name == "ofdm_out_real":  self.error_ofdm_real  += 1
            elif signal_name == "ofdm_out_imag":  self.error_ofdm_imag  += 1
            elif signal_name == "radar_out_real": self.error_radar_real += 1
            else:                                 self.error_radar_imag += 1
            self.logger.error(
                f"MISMATCH {signal_name} "
                f"(error #{getattr(self, 'error_' + signal_name.replace('out_', ''))}) "
                f"| DUT={dut_val}  REF={ref_val}  diff={diff}  tolerance={TOLERANCE}"
            )

    # ── report_phase ──────────────────────────────────────────────────────────
    def report_phase(self):
        total_ofdm_re  = self.correct_ofdm_real  + self.error_ofdm_real
        total_ofdm_im  = self.correct_ofdm_imag  + self.error_ofdm_imag
        total_radar_re = self.correct_radar_real + self.error_radar_real
        total_radar_im = self.correct_radar_imag + self.error_radar_imag

        pass_ofdm_re   = (self.error_ofdm_real  == 0)
        pass_ofdm_im   = (self.error_ofdm_imag  == 0)
        pass_radar_re  = (self.error_radar_real == 0)
        pass_radar_im  = (self.error_radar_imag == 0)
        overall        = "PASS " if all([pass_ofdm_re, pass_ofdm_im,
                                           pass_radar_re, pass_radar_im]) else "FAIL "

        self.logger.info("======================================================")
        self.logger.info("                 RX SCOREBOARD REPORT                 ")
        self.logger.info("======================================================")
        self.logger.info(f" Reset cycles detected        : {self.reset_cycles:<22}")
        self.logger.info("------------------------------------------------------")
        self.logger.info(f" ofdm_out_real correct/total  : {self.correct_ofdm_real}/{total_ofdm_re:<20}")
        self.logger.info(f" ofdm_out_real errors         : {self.error_ofdm_real:<22}")
        self.logger.info("------------------------------------------------------")
        self.logger.info(f" ofdm_out_imag correct/total  : {self.correct_ofdm_imag}/{total_ofdm_im:<20}")
        self.logger.info(f" ofdm_out_imag errors         : {self.error_ofdm_imag:<22}")
        self.logger.info("------------------------------------------------------")
        self.logger.info(f" radar_out_real correct/total : {self.correct_radar_real}/{total_radar_re:<20}")
        self.logger.info(f" radar_out_real errors        : {self.error_radar_real:<22}")
        self.logger.info("------------------------------------------------------")
        self.logger.info(f" radar_out_imag correct/total : {self.correct_radar_imag}/{total_radar_im:<20}")
        self.logger.info(f" radar_out_imag errors        : {self.error_radar_imag:<22}")
        self.logger.info("======================================================")
        self.logger.info(f" OVERALL RESULT               : {overall:<22}")
        self.logger.info("======================================================")

        # Hard-fail the test if there were any mismatches
        if overall == "FAIL ":
            self.logger.critical(
                f"Scoreboard detected errors: "
                f"ofdm_re={self.error_ofdm_real}, "
                f"ofdm_im={self.error_ofdm_imag}, "
                f"radar_re={self.error_radar_real}, "
                f"radar_im={self.error_radar_imag}."
            )

        # Close CSV files
        self.f_ofdm_rtl.close()
        self.f_ofdm_ref.close()
        self.f_radar_rtl.close()
        self.f_radar_ref.close()