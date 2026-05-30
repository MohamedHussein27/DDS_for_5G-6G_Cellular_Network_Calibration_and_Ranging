"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: system_top_scoreboard.py

    Description:
        Unified scoreboard for the full ISAC system (TX directly wired to RX).

        Mirrors the exact style of tx_scoreboard.py and rx_scoreboard.py.

        How it triggers the golden model
        ─────────────────────────────────
        Identical to the TX scoreboard: the system golden model needs only
        (FTW_start, FTW_step, cycles, ofdm_re_array, ofdm_im_array) to
        produce ALL output vectors in one call.  No sample accumulation is
        needed.  The scoreboard calls run_system_top_pipeline() on the first
        dds_item.enable pulse, then indexes into the three pre-computed
        reference vectors as each valid pulse arrives.

        Since TX is directly wired to RX with NO channel, ref_wr_en is never
        driven from outside.  The system golden model builds the reference RAM
        internally from the TX spectrum, so this scoreboard does NOT need to
        accumulate ref_wr_en writes.

        Three output paths checked
        ──────────────────────────
        1. TX output path
               DUT signal  : tx_valid + tx_out_re + tx_out_im
               Reference   : sys_golden_tx_re / sys_golden_tx_im   [4096]

        2. OFDM (communication) output path
               DUT signal  : ofdm_valid_out + ofdm_out_re + ofdm_out_im
               Reference   : sys_golden_ofdm_re / sys_golden_ofdm_im [2048]

        3. Radar (range profile) output path
               DUT signal  : radar_valid_out + radar_out_re + radar_out_im
               Reference   : sys_golden_radar_re / sys_golden_radar_im [2048]

        Tolerance: ±TOLERANCE LSBs on every comparison (same as TX and RX
        sub-scoreboards).

        CSV dump
        ────────
        Six CSV files are written during the simulation and closed in
        report_phase — one RTL and one Python reference per output path:
            tx_rtl_out.csv      /  tx_python_out.csv
            ofdm_rtl_out.csv    /  ofdm_python_out.csv
            radar_rtl_out.csv   /  radar_python_out.csv
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import logging
import numpy as np

# TX-side seq items (carry FTW params + tx_valid + ofdm_valid + radar_valid)
from system_top_seq_item import system_top_item
from dds_seq_item import *

# System golden model — single call produces all three output vectors
from system_top_golden_model import run_system_top_pipeline

# ─────────────────────────────────────────────────────────────────────────────
# RTL parameters
# ─────────────────────────────────────────────────────────────────────────────
WL          = 16
N_TX        = 4096    # TX output samples
N_RX        = 2048    # OFDM / Radar output samples each
TOLERANCE   = 4       # ±LSBs allowed between DUT and golden model


# ─────────────────────────────────────────────────────────────────────────────
# Fixed-point helpers  (identical style to both sub-scoreboards)
# ─────────────────────────────────────────────────────────────────────────────
def _sign16(v):
    """Sign-extend a raw 16-bit integer value from the DUT."""
    v = int(v) & 0xFFFF
    return v - 0x10000 if (v & 0x8000) else v


# ─────────────────────────────────────────────────────────────────────────────
# Scoreboard
# ─────────────────────────────────────────────────────────────────────────────
class system_top_scoreboard(uvm_scoreboard):
    """
    Unified scoreboard for the full ISAC system.
    TX is directly wired to RX — one golden model call covers everything.
    """

    def __init__(self, name, parent):
        super().__init__(name, parent)

    # ── build_phase ───────────────────────────────────────────────────────────
    def build_phase(self):
        # ── Analysis ports (same pattern as tx_scoreboard) ────────────────
        # Top-level system item (TX+RX DUT signals in one transaction)
        self.sb_export  = uvm_analysis_export("sb_export",  self)
        self.sb_fifo    = uvm_tlm_analysis_fifo("sb_fifo",  self)
        self.sb_export  = self.sb_fifo.analysis_export

        # DDS config item (FTW_start, FTW_step, cycles, enable)
        self.dds_export = uvm_analysis_export("dds_export", self)
        self.dds_fifo   = uvm_tlm_analysis_fifo("dds_fifo", self)
        self.dds_export = self.dds_fifo.analysis_export

        # DUT handle — used to backdoor-read the OFDM ROM (same as tx_scoreboard)
        self.dut = ConfigDB().get(self, "", "DUT")

        # ── Pass/fail counters — one pair per signal per path ─────────────
        # TX path
        self.correct_tx_real    = 0;  self.error_tx_real    = 0
        self.correct_tx_imag    = 0;  self.error_tx_imag    = 0
        # OFDM path
        self.correct_ofdm_real  = 0;  self.error_ofdm_real  = 0
        self.correct_ofdm_imag  = 0;  self.error_ofdm_imag  = 0
        # Radar path
        self.correct_radar_real = 0;  self.error_radar_real = 0
        self.correct_radar_imag = 0;  self.error_radar_imag = 0
        # Misc
        self.reset_cycles       = 0

        # ── Golden reference vectors (filled once on first dds_enable) ────
        self._ref_tx_real    = []
        self._ref_tx_imag    = []
        self._ref_ofdm_real  = []
        self._ref_ofdm_imag  = []
        self._ref_radar_real = []
        self._ref_radar_imag = []
        self._golden_ready   = False

        # ── Independent output-sample indices ─────────────────────────────
        self._tx_idx    = 0   # counts against N_TX  (4096)
        self._ofdm_idx  = 0   # counts against N_RX  (2048)
        self._radar_idx = 0   # counts against N_RX  (2048)

        # ── Global cycle counter (heartbeat / boundary prints) ────────────
        self._cycle = 0

        # ── CSV dump files ────────────────────────────────────────────────
        self.f_tx_rtl    = open("tx_rtl_out.csv",     "w")
        self.f_tx_ref    = open("tx_python_out.csv",   "w")
        self.f_ofdm_rtl  = open("ofdm_rtl_out.csv",   "w")
        self.f_ofdm_ref  = open("ofdm_python_out.csv", "w")
        self.f_radar_rtl = open("radar_rtl_out.csv",  "w")
        self.f_radar_ref = open("radar_python_out.csv","w")

        # Write CSV headers
        self.f_tx_rtl.write("RTL_Real,RTL_Imag\n")
        self.f_tx_ref.write("REF_Real,REF_Imag\n")
        self.f_ofdm_rtl.write("RTL_Real,RTL_Imag\n")
        self.f_ofdm_ref.write("REF_Real,REF_Imag\n")
        self.f_radar_rtl.write("RTL_Real,RTL_Imag\n")
        self.f_radar_ref.write("REF_Real,REF_Imag\n")

    # ── run_phase ─────────────────────────────────────────────────────────────
    async def run_phase(self):
        while True:
            item     = await self.sb_fifo.get()
            dds_item = await self.dds_fifo.get()
            self._process(item, dds_item)

    # ── _process ──────────────────────────────────────────────────────────────
    def _process(self, item, dds_item):

        # ── Reset: flush all state ────────────────────────────────────────
        if not item.rst_n:
            self.reset_cycles    += 1
            self._ref_tx_real     = []
            self._ref_tx_imag     = []
            self._ref_ofdm_real   = []
            self._ref_ofdm_imag   = []
            self._ref_radar_real  = []
            self._ref_radar_imag  = []
            self._golden_ready    = False
            self._tx_idx          = 0
            self._ofdm_idx        = 0
            self._radar_idx       = 0
            self._cycle           = 0
            self.logger.info("Reset detected — system golden model flushed.")
            return

        # ── Call system golden model once on first dds_enable ─────────────
        if dds_item.enable and not self._golden_ready:
            self.logger.info(
                f"Calling system golden model: "
                f"FTW_start={dds_item.FTW_start}, "
                f"FTW_step={dds_item.FTW_step}, "
                f"cycles={dds_item.cycles}"
            )

            # ── Backdoor-read OFDM ROM (identical to tx_scoreboard) ───────
            ofdm_re_list = [0] * N_RX
            ofdm_im_list = [0] * N_RX
            dut_ram_re   = self.dut.u_ofdm_rom.rom_real
            dut_ram_im   = self.dut.u_ofdm_rom.rom_imag
            ram_depth    = len(dut_ram_re)

            self.logger.info(
                f"Reading OFDM ROM (depth={ram_depth}) for system golden model..."
            )
            for mem_idx in range(ram_depth):
                try:
                    val_re = int(dut_ram_re[mem_idx].value)
                    val_im = int(dut_ram_im[mem_idx].value)
                except ValueError:
                    val_re = 0
                    val_im = 0
                ofdm_re_list[mem_idx] = _sign16(val_re)
                ofdm_im_list[mem_idx] = _sign16(val_im)

            # ── Single call: TX→RX system golden model ────────────────────
            # Returns all six output vectors at once.
            # ref_ram is built internally from the TX spectrum — no
            # ref_wr_en accumulation needed by this scoreboard.
            (tx_re_vec,    tx_im_vec,
             ofdm_re_vec,  ofdm_im_vec,
             radar_re_vec, radar_im_vec) = run_system_top_pipeline(
                FTW_start     = dds_item.FTW_start,
                FTW_step      = dds_item.FTW_step,
                N_cycles      = dds_item.cycles,
                Fs            = 491.52e6,
                ofdm_re_array = ofdm_re_list,
                ofdm_im_array = ofdm_im_list,
                debug_xlsx    = None,   # disable per-call dump; run manually if needed
            )

            self._ref_tx_real    = list(tx_re_vec)
            self._ref_tx_imag    = list(tx_im_vec)
            self._ref_ofdm_real  = list(ofdm_re_vec)
            self._ref_ofdm_imag  = list(ofdm_im_vec)
            self._ref_radar_real = list(radar_re_vec)
            self._ref_radar_imag = list(radar_im_vec)
            self._golden_ready   = True
            self._tx_idx         = 0
            self._ofdm_idx       = 0
            self._radar_idx      = 0

            # Boundary debug prints (mirrors tx_scoreboard style)
            print(f"Golden TX     sample #0:    re={self._ref_tx_real[0]}     im={self._ref_tx_imag[0]}")
            print(f"Golden TX     sample #4095: re={self._ref_tx_real[4095]}  im={self._ref_tx_imag[4095]}")
            print(f"Golden OFDM   sample #0:    re={self._ref_ofdm_real[0]}   im={self._ref_ofdm_imag[0]}")
            print(f"Golden OFDM   sample #2047: re={self._ref_ofdm_real[2047]} im={self._ref_ofdm_imag[2047]}")
            print(f"Golden Radar  sample #0:    re={self._ref_radar_real[0]}  im={self._ref_radar_imag[0]}")
            print(f"Golden Radar  sample #2047: re={self._ref_radar_real[2047]} im={self._ref_radar_imag[2047]}")

        # ── Heartbeat ─────────────────────────────────────────────────────
        if self._cycle % 1000 == 0:
            self.logger.info(
                f"--- SIMULATION HEARTBEAT: Processing cycle {self._cycle} ---"
            )

        # ── PATH 1: Compare TX output ─────────────────────────────────────
        if item.tx_valid:
            if not self._golden_ready:
                self.logger.warning(
                    "DUT asserted tx_valid but golden model not ready — skipping."
                )
            elif self._tx_idx >= N_TX:
                self.logger.warning(
                    f"tx_valid beyond frame size N={N_TX} "
                    f"(tx_idx={self._tx_idx}) — skipping."
                )
            else:
                dut_real = _sign16(item.tx_out_real)
                dut_imag = _sign16(item.tx_out_imag)
                ref_real = self._ref_tx_real[self._tx_idx]
                ref_imag = self._ref_tx_imag[self._tx_idx]

                if self._tx_idx < 4 or self._tx_idx > N_TX - 4:
                    self.logger.info(
                        f"TX compare sample #{self._tx_idx}: "
                        f"DUT(re={dut_real}, im={dut_imag}) vs "
                        f"REF(re={ref_real}, im={ref_imag})"
                    )

                self._compare("tx_out_real", dut_real, ref_real)
                self._compare("tx_out_imag", dut_imag, ref_imag)

                self.f_tx_rtl.write(f"{dut_real},{dut_imag}\n")
                self.f_tx_ref.write(f"{ref_real},{ref_imag}\n")
                self._tx_idx += 1

        # ── PATH 2: Compare OFDM output ───────────────────────────────────
        if item.ofdm_valid_out:
            if not self._golden_ready:
                self.logger.warning(
                    "DUT asserted ofdm_valid_out but golden model not ready — skipping."
                )
            elif self._ofdm_idx >= N_RX:
                self.logger.warning(
                    f"ofdm_valid_out beyond frame size N={N_RX} "
                    f"(ofdm_idx={self._ofdm_idx}) — skipping."
                )
            else:
                dut_real = _sign16(item.ofdm_out_re)
                dut_imag = _sign16(item.ofdm_out_im)
                ref_real = self._ref_ofdm_real[self._ofdm_idx]
                ref_imag = self._ref_ofdm_imag[self._ofdm_idx]

                if self._ofdm_idx < 4 or self._ofdm_idx > N_RX - 4:
                    self.logger.info(
                        f"OFDM compare sample #{self._ofdm_idx}: "
                        f"DUT(re={dut_real}, im={dut_imag}) vs "
                        f"REF(re={ref_real}, im={ref_imag})"
                    )

                self._compare("ofdm_out_real", dut_real, ref_real)
                self._compare("ofdm_out_imag", dut_imag, ref_imag)

                self.f_ofdm_rtl.write(f"{dut_real},{dut_imag}\n")
                self.f_ofdm_ref.write(f"{ref_real},{ref_imag}\n")
                self._ofdm_idx += 1

        # ── PATH 3: Compare Radar output ──────────────────────────────────
        if item.radar_valid_out:
            if not self._golden_ready:
                self.logger.warning(
                    "DUT asserted radar_valid_out but golden model not ready — skipping."
                )
            elif self._radar_idx >= N_RX:
                self.logger.warning(
                    f"radar_valid_out beyond frame size N={N_RX} "
                    f"(radar_idx={self._radar_idx}) — skipping."
                )
            else:
                dut_real = _sign16(item.radar_out_re)
                dut_imag = _sign16(item.radar_out_im)
                ref_real = self._ref_radar_real[self._radar_idx]
                ref_imag = self._ref_radar_imag[self._radar_idx]

                if self._radar_idx < 4 or self._radar_idx > N_RX - 4:
                    self.logger.info(
                        f"Radar compare sample #{self._radar_idx}: "
                        f"DUT(re={dut_real}, im={dut_imag}) vs "
                        f"REF(re={ref_real}, im={ref_imag})"
                    )

                self._compare("radar_out_real", dut_real, ref_real)
                self._compare("radar_out_imag", dut_imag, ref_imag)

                self.f_radar_rtl.write(f"{dut_real},{dut_imag}\n")
                self.f_radar_ref.write(f"{ref_real},{ref_imag}\n")
                self._radar_idx += 1

        self._cycle += 1

    # ── _compare ──────────────────────────────────────────────────────────────
    def _compare(self, signal_name: str, dut_val: int, ref_val: int):
        """Compare one signal with ±TOLERANCE LSB window."""
        diff = abs(dut_val - ref_val)
        if diff <= TOLERANCE:
            if   signal_name == "tx_out_real":    self.correct_tx_real    += 1
            elif signal_name == "tx_out_imag":    self.correct_tx_imag    += 1
            elif signal_name == "ofdm_out_real":  self.correct_ofdm_real  += 1
            elif signal_name == "ofdm_out_imag":  self.correct_ofdm_imag  += 1
            elif signal_name == "radar_out_real": self.correct_radar_real += 1
            else:                                 self.correct_radar_imag += 1
        else:
            if   signal_name == "tx_out_real":    self.error_tx_real    += 1
            elif signal_name == "tx_out_imag":    self.error_tx_imag    += 1
            elif signal_name == "ofdm_out_real":  self.error_ofdm_real  += 1
            elif signal_name == "ofdm_out_imag":  self.error_ofdm_imag  += 1
            elif signal_name == "radar_out_real": self.error_radar_real += 1
            else:                                 self.error_radar_imag += 1
            self.logger.error(
                f"MISMATCH {signal_name} "
                f"| DUT={dut_val}  REF={ref_val}  diff={diff}  tol={TOLERANCE}"
            )

    # ── report_phase ──────────────────────────────────────────────────────────
    def report_phase(self):
        # ── Totals ────────────────────────────────────────────────────────
        tot_tx_re    = self.correct_tx_real    + self.error_tx_real
        tot_tx_im    = self.correct_tx_imag    + self.error_tx_imag
        tot_ofdm_re  = self.correct_ofdm_real  + self.error_ofdm_real
        tot_ofdm_im  = self.correct_ofdm_imag  + self.error_ofdm_imag
        tot_radar_re = self.correct_radar_real + self.error_radar_real
        tot_radar_im = self.correct_radar_imag + self.error_radar_imag

        # ── Per-path pass/fail ────────────────────────────────────────────
        pass_tx    = (self.error_tx_real    == 0 and self.error_tx_imag    == 0)
        pass_ofdm  = (self.error_ofdm_real  == 0 and self.error_ofdm_imag  == 0)
        pass_radar = (self.error_radar_real == 0 and self.error_radar_imag == 0)
        overall    = "PASS ✓" if (pass_tx and pass_ofdm and pass_radar) else "FAIL ✗"

        def _pf(b): return "PASS ✓" if b else "FAIL ✗"

        self.logger.info("##############################################################")
        self.logger.info("           SYSTEM TOP SCOREBOARD REPORT               ")
        self.logger.info("##############################################################")
        self.logger.info(f"  Reset cycles detected          : {self.reset_cycles:<20} ")
        self.logger.info("###################### TX PATH #################################")
        self.logger.info(f"  tx_out_real  correct / total   : {self.correct_tx_real}/{tot_tx_re:<18} ")
        self.logger.info(f"  tx_out_real  errors            : {self.error_tx_real:<20} ")
        self.logger.info(f"  tx_out_imag  correct / total   : {self.correct_tx_imag}/{tot_tx_im:<18} ")
        self.logger.info(f"  tx_out_imag  errors            : {self.error_tx_imag:<20} ")
        self.logger.info(f"  TX PATH RESULT                 : {_pf(pass_tx):<20} ")
        self.logger.info("############################### OFDM PATH ###############################")
        self.logger.info(f"  ofdm_out_real  correct / total : {self.correct_ofdm_real}/{tot_ofdm_re:<18} ")
        self.logger.info(f"  ofdm_out_real  errors          : {self.error_ofdm_real:<20} ")
        self.logger.info(f"  ofdm_out_imag  correct / total : {self.correct_ofdm_imag}/{tot_ofdm_im:<18} ")
        self.logger.info(f"  ofdm_out_imag  errors          : {self.error_ofdm_imag:<20} ")
        self.logger.info(f"  OFDM PATH RESULT               : {_pf(pass_ofdm):<20} ")
        self.logger.info("############################### RADAR PATH ###############################")
        self.logger.info(f"  radar_out_real correct / total : {self.correct_radar_real}/{tot_radar_re:<18} ")
        self.logger.info(f"  radar_out_real errors          : {self.error_radar_real:<20} ")
        self.logger.info(f"  radar_out_imag correct / total : {self.correct_radar_imag}/{tot_radar_im:<18} ")
        self.logger.info(f"  radar_out_imag errors          : {self.error_radar_imag:<20} ")
        self.logger.info(f"  RADAR PATH RESULT              : {_pf(pass_radar):<20} ")
        self.logger.info("###############################")
        self.logger.info(f"  OVERALL RESULT                 : {overall:<20} ")
        self.logger.info("###############################")

        # Hard-fail if any path has errors
        if overall == "FAIL ✗":
            self.logger.critical(
                f"Scoreboard errors — "
                f"TX: re={self.error_tx_real} im={self.error_tx_imag} | "
                f"OFDM: re={self.error_ofdm_real} im={self.error_ofdm_imag} | "
                f"Radar: re={self.error_radar_real} im={self.error_radar_imag}"
            )

        # ── Close CSV files ───────────────────────────────────────────────
        self.f_tx_rtl.close();    self.f_tx_ref.close()
        self.f_ofdm_rtl.close();  self.f_ofdm_ref.close()
        self.f_radar_rtl.close(); self.f_radar_ref.close()