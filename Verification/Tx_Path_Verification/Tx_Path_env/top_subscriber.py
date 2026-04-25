"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: top_subscriber.py

    Description:
        This functional coverage subscriber passively monitors the TX_TOP boundary 
        to ensure the system-level verification environment adequately stimulates 
        the top-level datapath and control signals. It tracks DDS controls, OFDM 
        data path handshaking (ofdm_rd_en, tx_valid), and ensures the 16-bit 
        signed fixed-point data extensively hits critical boundary conditions 
        for both the input OFDM signals and the final TX output.
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from cocotb_coverage.coverage import CoverPoint, CoverCross, coverage_section
from cocotb.queue import *
from cocotb_coverage.coverage import coverage_db
from pyuvm import *
import pyuvm
import logging
 
from seq_item import top_item
 
# ─────────────────────────────────────────────────────────────────────────────
# Constants that mirror the RTL parameters
# ─────────────────────────────────────────────────────────────────────────────
WL   = 16                         # System Word Length (signed)
N    = 4096                       # FFT/IFFT size
MAX_VAL  =  (1 << (WL - 1)) - 1   #  32767
MIN_VAL  = -(1 << (WL - 1))       # -32768
MID_POS  =  MAX_VAL // 2          #  16383
MID_NEG  =  MIN_VAL // 2          # -16384
 
# ─────────────────────────────────────────────────────────────────────────────
# Helper – sign-extend a raw integer read from a 16-bit RTL signal
# ─────────────────────────────────────────────────────────────────────────────
def _sign16(v):
    v = int(v) & 0xFFFF
    if v & 0x8000:
        return v - 0x10000
    else:
        return v
 
 
class top_subscriber(uvm_subscriber):
    """
    Functional-coverage subscriber for the TX_TOP system.
 
    Cover points
    ─────────────
    CP1  dds_enable      : both de-asserted (0) and asserted (1)
    CP2  ofdm_rd_en      : both de-asserted (0) and asserted (1)
    CP3  tx_valid        : both de-asserted (0) and asserted (1)
    CP4  rst_n           : reset active (0) and inactive (1)
    CP5  ofdm_in_real    : negative / zero / positive / max-positive / min-negative
    CP6  ofdm_in_imag    : negative / zero / positive / max-positive / min-negative
    CP7  tx_out_real     : negative / zero / positive / max-positive / min-negative
    CP8  tx_out_imag     : negative / zero / positive / max-positive / min-negative
    CP9  ofdm_in quadrant cross : (pos,pos) (pos,neg) (neg,pos) (neg,neg)
    CP10 tx_out quadrant cross  : (pos,pos) (pos,neg) (neg,pos) (neg,neg) cross with tx_valid
    """
 
    def __init__(self, name, parent):
        super().__init__(name, parent)
 
    # ──────────────────────────────────────────────────────────────────────
    # build_phase : create export + FIFO 
    # ──────────────────────────────────────────────────────────────────────
    def build_phase(self):
        self.sub_export = uvm_analysis_export("sub_export", self)
        self.sub_fifo   = uvm_tlm_analysis_fifo("sub_fifo", self)
        # re-point the export straight at the FIFO's analysis_export
        self.sub_export = self.sub_fifo.analysis_export # back to this 
 
    # ──────────────────────────────────────────────────────────────────────
    # write : mandatory uvm_subscriber hook
    # ──────────────────────────────────────────────────────────────────────
    def write(self, item):
        self.logger.debug(f"Top Subscriber received: {item.convert2string()}")
 
    # ══════════════════════════════════════════════════════════════════════
    # COVER POINTS  (all bins defined with lambda pin functions)
    # ══════════════════════════════════════════════════════════════════════
 
    # ── CP1 : dds_enable ──────────────────────────────────────────────────
    @CoverPoint(
        "top.dds_enable",
        xf   = lambda tr: int(tr.dds_enable),
        bins = [0, 1],
        bins_labels = ["dds_disabled", "dds_enabled"],
    )
    # ── CP2 : ofdm_rd_en ──────────────────────────────────────────────────
    @CoverPoint(
        "top.ofdm_rd_en",
        xf   = lambda tr: int(tr.ofdm_rd_en),
        bins = [0, 1],
        bins_labels = ["rd_en_deasserted", "rd_en_asserted"],
    )
    # ── CP3 : tx_valid ────────────────────────────────────────────────────
    @CoverPoint(
        "top.tx_valid",
        xf   = lambda tr: int(tr.tx_valid),
        bins = [0, 1],
        bins_labels = ["tx_valid_deasserted", "tx_valid_asserted"],
    )
    # ── CP4 : rst_n ───────────────────────────────────────────────────────
    @CoverPoint(
        "top.rst_n",
        xf   = lambda tr: int(tr.rst_n),
        bins = [0, 1],
        bins_labels = ["reset_active", "reset_inactive"],
    )
    # ── CP5 : ofdm_in_real range ──────────────────────────────────────────
    @CoverPoint(
        "top.ofdm_in_real_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.ofdm_in_real) == MIN_VAL  else
            "negative" if _sign16(tr.ofdm_in_real) <  0        else
            "zero"     if _sign16(tr.ofdm_in_real) == 0        else
            "max_pos"  if _sign16(tr.ofdm_in_real) == MAX_VAL  else
            "positive"
        ),
        bins        = ["min_neg", "negative", "zero", "positive", "max_pos"],
        bins_labels = ["ofdm_re_min_neg", "ofdm_re_neg",
                       "ofdm_re_zero", "ofdm_re_pos", "ofdm_re_max_pos"],
    )
    # ── CP6 : ofdm_in_imag range ──────────────────────────────────────────
    @CoverPoint(
        "top.ofdm_in_imag_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.ofdm_in_imag) == MIN_VAL  else
            "negative" if _sign16(tr.ofdm_in_imag) <  0        else
            "zero"     if _sign16(tr.ofdm_in_imag) == 0        else
            "max_pos"  if _sign16(tr.ofdm_in_imag) == MAX_VAL  else
            "positive"
        ),
        bins        = ["min_neg", "negative", "zero", "positive", "max_pos"],
        bins_labels = ["ofdm_im_min_neg", "ofdm_im_neg",
                       "ofdm_im_zero", "ofdm_im_pos", "ofdm_im_max_pos"],
    )
    # ── CP7 : tx_out_real range ───────────────────────────────────────────
    @CoverPoint(
        "top.tx_out_real_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.tx_out_real) == MIN_VAL  else
            "negative" if _sign16(tr.tx_out_real) <  0        else
            "zero"     if _sign16(tr.tx_out_real) == 0        else
            "max_pos"  if _sign16(tr.tx_out_real) == MAX_VAL  else
            "positive"
        ),
        bins        = ["min_neg", "negative", "zero", "positive", "max_pos"],
        bins_labels = ["tx_re_min_neg", "tx_re_neg",
                       "tx_re_zero", "tx_re_pos", "tx_re_max_pos"],
    )
    # ── CP8 : tx_out_imag range ───────────────────────────────────────────
    @CoverPoint(
        "top.tx_out_imag_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.tx_out_imag) == MIN_VAL  else
            "negative" if _sign16(tr.tx_out_imag) <  0        else
            "zero"     if _sign16(tr.tx_out_imag) == 0        else
            "max_pos"  if _sign16(tr.tx_out_imag) == MAX_VAL  else
            "positive"
        ),
        bins        = ["min_neg", "negative", "zero", "positive", "max_pos"],
        bins_labels = ["tx_im_min_neg", "tx_im_neg",
                       "tx_im_zero", "tx_im_pos", "tx_im_max_pos"],
    )
    # ── CP9 : ofdm_in_real × ofdm_in_imag quadrant cross ──────────────────
    @CoverCross(
        "top.ofdm_input_quadrant",
        items = ["top.ofdm_in_real_range", "top.ofdm_in_imag_range"],
    )
    # ── CP10 : tx_out_real × tx_out_imag quadrant, tx_valid cross ─────────
    @CoverCross(
        "top.tx_output_quadrant_valid",
        items = ["top.tx_valid", "top.tx_out_real_range", "top.tx_out_imag_range"],
    )

    def sample(self, tr):
        pass
 
    # ──────────────────────────────────────────────────────────────────────
    # run_phase : drain the FIFO and sample every item
    # ──────────────────────────────────────────────────────────────────────
    async def run_phase(self):
        while True:
            seq_item_sub = await self.sub_fifo.get()
            self.sample(seq_item_sub)