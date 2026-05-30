"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: system_top_subscriber.py

    Description:
        Functional-coverage subscriber for the full ISAC SYSTEM_TOP.

        Receives system_top_item transactions from system_top_monitor
        and samples them against all cover points / crosses.

        The covergroup mirrors and extends both tx_top_subscriber.py and
        rx_top_subscriber.py, adding system-level crosses that can only
        be verified at this boundary (e.g. simultaneous TX and RX valid).

        Cover points
        ─────────────
        TX Control
          CP01  rst_n             : active (0) / inactive (1)
          CP02  wr_en             : 0 / 1
          CP03  rd_en             : 0 / 1
          CP04  addr              : 0x0 (FTW_start) / 0x4 (FTW_step) / 0x8 (cycles)
          CP05  dds_ready_flag    : 0 / 1
          CP06  tx_valid          : 0 / 1

        TX Output Data Ranges (16-bit signed: min_neg/negative/zero/positive/max_pos)
          CP07  tx_out_real_range
          CP08  tx_out_imag_range

        RX Loopback Input
          CP09  rx_valid_in       : 0 / 1

        OFDM Output
          CP10  ofdm_valid_out    : 0 / 1
          CP11  ofdm_out_real_range
          CP12  ofdm_out_imag_range

        Radar Output
          CP13  radar_valid_out   : 0 / 1
          CP14  radar_out_real_range
          CP15  radar_out_imag_range

        Crosses
          CP16  bus_write_cross   : addr × wr_en  (ensure every address was written)
          CP17  bus_read_cross    : addr × rd_en  (ensure every address was read)
          CP18  tx_output_quadrant_valid  : tx_valid × tx_re_range × tx_im_range
          CP19  ofdm_output_quadrant_valid: ofdm_valid × ofdm_re × ofdm_im
          CP20  radar_output_quadrant_valid: radar_valid × radar_re × radar_im
          CP21  system_concurrent_valid  : tx_valid × ofdm_valid × radar_valid
                (verifies all three output paths fire simultaneously)
"""

import cocotb
from cocotb.triggers import *
from cocotb_coverage.crv import *
from cocotb_coverage.coverage import CoverPoint, CoverCross
from cocotb_coverage.coverage import coverage_db
from pyuvm import *
import pyuvm
import logging

from system_top_seq_item import system_top_item

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
WL      = 16
MAX_VAL =  (1 << (WL - 1)) - 1   #  32767
MIN_VAL = -(1 << (WL - 1))       # -32768


def _sign16(v):
    v = int(v) & 0xFFFF
    return v - 0x10000 if (v & 0x8000) else v


def _range_bucket(raw):
    """Map a raw 16-bit signed integer to a named range bucket."""
    v = _sign16(raw)
    if   v == MIN_VAL: return "min_neg"
    elif v <  0:       return "negative"
    elif v == 0:       return "zero"
    elif v == MAX_VAL: return "max_pos"
    else:              return "positive"


_RANGE_BINS = ["min_neg", "negative", "zero", "positive", "max_pos"]


class system_top_subscriber(uvm_subscriber):

    def __init__(self, name, parent):
        super().__init__(name, parent)

    # ── build_phase ───────────────────────────────────────────────────
    def build_phase(self):
        super().build_phase()
        self.sub_export = uvm_analysis_export("sub_export", self)
        self.sub_fifo   = uvm_tlm_analysis_fifo("sub_fifo", self)
        # Re-point the export at the FIFO's analysis_export
        self.sub_export = self.sub_fifo.analysis_export

    # ── write ─────────────────────────────────────────────────────────
    def write(self, item):
        """Called instantly by the monitor's analysis port."""
        # Skip pure-software backdoor transactions
        if hasattr(item, "is_backdoor") and item.is_backdoor:
            return
        self.sample(item)

    # 
    # COVER POINTS
    # 

    # ── CP01 : rst_n ──────────────────────────────────────────────────
    @CoverPoint(
        "sys.rst_n",
        xf   = lambda tr: int(tr.rst_n),
        bins = [0, 1],
        bins_labels = ["reset_active", "reset_inactive"],
    )
    # ── CP02 : wr_en ──────────────────────────────────────────────────
    @CoverPoint(
        "sys.wr_en",
        xf   = lambda tr: int(tr.wr_en),
        bins = [0, 1],
        bins_labels = ["wr_inactive", "wr_active"],
    )
    # ── CP03 : rd_en ──────────────────────────────────────────────────
    @CoverPoint(
        "sys.rd_en",
        xf   = lambda tr: int(tr.rd_en),
        bins = [0, 1],
        bins_labels = ["rd_inactive", "rd_active"],
    )
    # ── CP04 : addr ───────────────────────────────────────────────────
    @CoverPoint(
        "sys.addr",
        xf   = lambda tr: int(tr.addr),
        bins = [0x0, 0x4, 0x8],
        bins_labels = ["ADDR_FTW_START", "ADDR_FTW_STEP", "ADDR_CYCLES"],
    )
    # ── CP05 : dds_ready_flag ─────────────────────────────────────────
    @CoverPoint(
        "sys.dds_ready_flag",
        xf   = lambda tr: int(tr.dds_ready_flag),
        bins = [0, 1],
        bins_labels = ["dds_not_ready", "dds_ready"],
    )
    # ── CP06 : tx_valid ───────────────────────────────────────────────
    @CoverPoint(
        "sys.tx_valid",
        xf   = lambda tr: int(tr.tx_valid),
        bins = [0, 1],
        bins_labels = ["tx_valid_inactive", "tx_valid_active"],
    )
    # ── CP07 : tx_out_real range ──────────────────────────────────────
    @CoverPoint(
        "sys.tx_out_real_range",
        xf   = lambda tr: _range_bucket(tr.tx_out_real),
        bins = _RANGE_BINS,
    )
    # ── CP08 : tx_out_imag range ──────────────────────────────────────
    @CoverPoint(
        "sys.tx_out_imag_range",
        xf   = lambda tr: _range_bucket(tr.tx_out_imag),
        bins = _RANGE_BINS,
    )
    # ── CP09 : rx_valid_in ────────────────────────────────────────────
    @CoverPoint(
        "sys.rx_valid_in",
        xf   = lambda tr: int(tr.rx_valid_in),
        bins = [0, 1],
        bins_labels = ["rx_vld_inactive", "rx_vld_active"],
    )
    # ── CP10 : ofdm_valid_out ─────────────────────────────────────────
    @CoverPoint(
        "sys.ofdm_valid_out",
        xf   = lambda tr: int(tr.ofdm_valid_out),
        bins = [0, 1],
        bins_labels = ["ofdm_vld_inactive", "ofdm_vld_active"],
    )
    # ── CP11 : ofdm_out_real range ────────────────────────────────────
    @CoverPoint(
        "sys.ofdm_out_real_range",
        xf   = lambda tr: _range_bucket(tr.ofdm_out_re),
        bins = _RANGE_BINS,
    )
    # ── CP12 : ofdm_out_imag range ────────────────────────────────────
    @CoverPoint(
        "sys.ofdm_out_imag_range",
        xf   = lambda tr: _range_bucket(tr.ofdm_out_im),
        bins = _RANGE_BINS,
    )
    # ── CP13 : radar_valid_out ────────────────────────────────────────
    @CoverPoint(
        "sys.radar_valid_out",
        xf   = lambda tr: int(tr.radar_valid_out),
        bins = [0, 1],
        bins_labels = ["radar_vld_inactive", "radar_vld_active"],
    )
    # ── CP14 : radar_out_real range ───────────────────────────────────
    @CoverPoint(
        "sys.radar_out_real_range",
        xf   = lambda tr: _range_bucket(tr.radar_out_re),
        bins = _RANGE_BINS,
    )
    # ── CP15 : radar_out_imag range ───────────────────────────────────
    @CoverPoint(
        "sys.radar_out_imag_range",
        xf   = lambda tr: _range_bucket(tr.radar_out_im),
        bins = _RANGE_BINS,
    )

    # 
    # CROSSES
    # 

    # ── CP16 : ensure every address was written ───────────────────────
    @CoverCross(
        "sys.bus_write_cross",
        items    = ["sys.addr", "sys.wr_en"],
        ign_bins = [
            ("ADDR_FTW_START", "wr_inactive"),
            ("ADDR_FTW_STEP",  "wr_inactive"),
            ("ADDR_CYCLES",    "wr_inactive"),
        ],
    )
    # ── CP17 : ensure every address was read ──────────────────────────
    @CoverCross(
        "sys.bus_read_cross",
        items    = ["sys.addr", "sys.rd_en"],
        ign_bins = [
            ("ADDR_FTW_START", "rd_inactive"),
            ("ADDR_FTW_STEP",  "rd_inactive"),
            ("ADDR_CYCLES",    "rd_inactive"),
        ],
    )
    # ── CP18 : TX output quadrant × valid ─────────────────────────────
    @CoverCross(
        "sys.tx_output_quadrant_valid",
        items = ["sys.tx_valid", "sys.tx_out_real_range", "sys.tx_out_imag_range"],
    )
    # ── CP19 : OFDM output quadrant × valid ───────────────────────────
    @CoverCross(
        "sys.ofdm_output_quadrant_valid",
        items = ["sys.ofdm_valid_out", "sys.ofdm_out_real_range", "sys.ofdm_out_imag_range"],
    )
    # ── CP20 : Radar output quadrant × valid ──────────────────────────
    @CoverCross(
        "sys.radar_output_quadrant_valid",
        items = ["sys.radar_valid_out", "sys.radar_out_real_range", "sys.radar_out_imag_range"],
    )
    # ── CP21 : All three valid signals active simultaneously ──────────
    @CoverCross(
        "sys.system_concurrent_valid",
        items = ["sys.tx_valid", "sys.ofdm_valid_out", "sys.radar_valid_out"],
    )

    def sample(self, tr):
        pass   # triggered by decorators above

    # ── run_phase ─────────────────────────────────────────────────────
    async def run_phase(self):
        while True:
            item = await self.sub_fifo.get()
            self.sample(item)

    # ── report_phase ──────────────────────────────────────────────────
    def report_phase(self):
        super().report_phase()
        coverage_db.export_to_xml("system_top_coverage.xml")

        self.logger.info("###############################################################")
        self.logger.info("           SYSTEM TOP FUNCTIONAL COVERAGE REPORT              ")
        self.logger.info("################################################################")

        all_pass = True
        for name in coverage_db:
            if name.startswith("sys."):
                cp         = coverage_db[name]
                percentage = cp.cover_percentage
                status     = "PASS ✓" if percentage >= 100.0 else "FAIL ✗"
                if percentage < 100.0:
                    all_pass = False
                short = name.split(".")[-1]
                self.logger.info(
                    f"║  {short:<34} : {status:<8} ({percentage:>5.1f}%)  ║"
                )

        overall = "PASS " if all_pass else "FAIL "
        self.logger.info("#########################################################")
        self.logger.info(f"  OVERALL COVERAGE TARGETS           : {overall:<17}  ")
        self.logger.info("#########################################################")
        self.logger.info("Detailed coverage data exported to system_top_coverage.xml")
