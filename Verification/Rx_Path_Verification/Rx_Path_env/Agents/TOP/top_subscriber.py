"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: rx_subscriber.py

    Description:
        This functional coverage subscriber passively monitors the RX_TOP boundary 
        to ensure the system-level verification environment adequately stimulates 
        both the incoming channel data and the reference RAM. It tracks validation
        signals for the two distinct outputs (ofdm_out and radar_out) and ensures 
        the 16-bit signed fixed-point data extensively hits critical boundary 
        conditions across both processing paths.
"""
import cocotb
from cocotb.triggers import *
from cocotb_coverage.crv import *
from cocotb_coverage.coverage import CoverPoint, CoverCross, coverage_section
from cocotb_coverage.coverage import coverage_db
from pyuvm import *
import pyuvm
import logging
 
from top_seq_item import rx_item
 
# ─────────────────────────────────────────────────────────────────────────────
# Constants that mirror the RTL parameters
# ─────────────────────────────────────────────────────────────────────────────
WL   = 16                         # System Word Length (signed)
MAX_VAL  =  (1 << (WL - 1)) - 1   #  32767
MIN_VAL  = -(1 << (WL - 1))       # -32768
 
# ─────────────────────────────────────────────────────────────────────────────
# Helper – sign-extend a raw integer read from a 16-bit RTL signal
# ─────────────────────────────────────────────────────────────────────────────
def _sign16(v):
    v = int(v) & 0xFFFF
    if v & 0x8000:
        return v - 0x10000
    else:
        return v
 
 
class rx_subscriber(uvm_subscriber):
    """
    Functional-coverage subscriber for the dual-path RX_TOP system.
 
    Cover points
    ─────────────
    CP1  rx_valid_in      : both de-asserted (0) and asserted (1)
    CP2  ref_wr_en        : both de-asserted (0) and asserted (1)
    CP3  ofdm_valid_out   : both de-asserted (0) and asserted (1)
    CP4  radar_valid_out  : both de-asserted (0) and asserted (1)
    CP5  rst_n            : reset active (0) and inactive (1)
    
    Data Boundary Coverpoints (Negative / Zero / Positive / Min-Neg / Max-Pos):
    CP6-7   rx_in (Real/Imag)
    CP8-9   ref_wr (Real/Imag)
    CP10-11 ofdm_out (Real/Imag)
    CP12-13 radar_out (Real/Imag)
    
    Crosses:
    CP14 rx_in quadrant
    CP15 ofdm_out quadrant crossed with ofdm_valid_out
    CP16 radar_out quadrant crossed with radar_valid_out
    """
 
    def __init__(self, name, parent):
        super().__init__(name, parent)
 
    # ──────────────────────────────────────────────────────────────────────
    # build_phase : create export + FIFO 
    # ──────────────────────────────────────────────────────────────────────
    def build_phase(self):
        super().build_phase()
        self.sub_export = uvm_analysis_export("sub_export", self)
        self.sub_fifo   = uvm_tlm_analysis_fifo("sub_fifo", self)
        # re-point the export straight at the FIFO's analysis_export
        self.sub_export = self.sub_fifo.analysis_export 
 
    # ──────────────────────────────────────────────────────────────────────
    # write : mandatory uvm_subscriber hook
    # ──────────────────────────────────────────────────────────────────────
    def write(self, item):
        self.logger.debug(f"RX Subscriber received: {item.convert2string()}")
 
    # ══════════════════════════════════════════════════════════════════════
    # COVER POINTS  (all bins defined with lambda pin functions)
    # ══════════════════════════════════════════════════════════════════════
 
    # ── Control Signals ───────────────────────────────────────────────────
    @CoverPoint("rx.rx_valid_in", xf = lambda tr: int(tr.rx_valid_in), bins = [0, 1])
    @CoverPoint("rx.ref_wr_en", xf = lambda tr: int(tr.ref_wr_en), bins = [0, 1])
    @CoverPoint("rx.ofdm_valid_out", xf = lambda tr: int(tr.ofdm_valid_out), bins = [0, 1])
    @CoverPoint("rx.radar_valid_out", xf = lambda tr: int(tr.radar_valid_out), bins = [0, 1])
    @CoverPoint("rx.rst_n", xf = lambda tr: int(tr.rst_n), bins = [0, 1])

    # ── Channel Input (rx_in) Data Ranges ─────────────────────────────────
    @CoverPoint(
        "rx.rx_in_real_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.rx_in_re) == MIN_VAL  else
            "negative" if _sign16(tr.rx_in_re) <  0        else
            "zero"     if _sign16(tr.rx_in_re) == 0        else
            "max_pos"  if _sign16(tr.rx_in_re) == MAX_VAL  else
            "positive"
        ),
        bins = ["min_neg", "negative", "zero", "positive", "max_pos"],
    )
    @CoverPoint(
        "rx.rx_in_imag_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.rx_in_im) == MIN_VAL  else
            "negative" if _sign16(tr.rx_in_im) <  0        else
            "zero"     if _sign16(tr.rx_in_im) == 0        else
            "max_pos"  if _sign16(tr.rx_in_im) == MAX_VAL  else
            "positive"
        ),
        bins = ["min_neg", "negative", "zero", "positive", "max_pos"],
    )

    # ── Reference RAM Input (ref_wr) Data Ranges ──────────────────────────
    @CoverPoint(
        "rx.ref_wr_real_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.ref_wr_re) == MIN_VAL  else
            "negative" if _sign16(tr.ref_wr_re) <  0        else
            "zero"     if _sign16(tr.ref_wr_re) == 0        else
            "max_pos"  if _sign16(tr.ref_wr_re) == MAX_VAL  else
            "positive"
        ),
        bins = ["min_neg", "negative", "zero", "positive", "max_pos"],
    )
    @CoverPoint(
        "rx.ref_wr_imag_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.ref_wr_im) == MIN_VAL  else
            "negative" if _sign16(tr.ref_wr_im) <  0        else
            "zero"     if _sign16(tr.ref_wr_im) == 0        else
            "max_pos"  if _sign16(tr.ref_wr_im) == MAX_VAL  else
            "positive"
        ),
        bins = ["min_neg", "negative", "zero", "positive", "max_pos"],
    )

    # ── OFDM Output (ofdm_out) Data Ranges ────────────────────────────────
    @CoverPoint(
        "rx.ofdm_out_real_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.ofdm_out_re) == MIN_VAL  else
            "negative" if _sign16(tr.ofdm_out_re) <  0        else
            "zero"     if _sign16(tr.ofdm_out_re) == 0        else
            "max_pos"  if _sign16(tr.ofdm_out_re) == MAX_VAL  else
            "positive"
        ),
        bins = ["min_neg", "negative", "zero", "positive", "max_pos"],
    )
    @CoverPoint(
        "rx.ofdm_out_imag_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.ofdm_out_im) == MIN_VAL  else
            "negative" if _sign16(tr.ofdm_out_im) <  0        else
            "zero"     if _sign16(tr.ofdm_out_im) == 0        else
            "max_pos"  if _sign16(tr.ofdm_out_im) == MAX_VAL  else
            "positive"
        ),
        bins = ["min_neg", "negative", "zero", "positive", "max_pos"],
    )

    # ── Radar Output (radar_out) Data Ranges ──────────────────────────────
    @CoverPoint(
        "rx.radar_out_real_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.radar_out_re) == MIN_VAL  else
            "negative" if _sign16(tr.radar_out_re) <  0        else
            "zero"     if _sign16(tr.radar_out_re) == 0        else
            "max_pos"  if _sign16(tr.radar_out_re) == MAX_VAL  else
            "positive"
        ),
        bins = ["min_neg", "negative", "zero", "positive", "max_pos"],
    )
    @CoverPoint(
        "rx.radar_out_imag_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.radar_out_im) == MIN_VAL  else
            "negative" if _sign16(tr.radar_out_im) <  0        else
            "zero"     if _sign16(tr.radar_out_im) == 0        else
            "max_pos"  if _sign16(tr.radar_out_im) == MAX_VAL  else
            "positive"
        ),
        bins = ["min_neg", "negative", "zero", "positive", "max_pos"],
    )

    # ── CROSS COVERAGE ────────────────────────────────────────────────────
    @CoverCross(
        "rx.channel_input_quadrant",
        items = ["rx.rx_in_real_range", "rx.rx_in_imag_range"],
    )
    @CoverCross(
        "rx.ofdm_output_quadrant_valid",
        items = ["rx.ofdm_valid_out", "rx.ofdm_out_real_range", "rx.ofdm_out_imag_range"],
    )
    @CoverCross(
        "rx.radar_output_quadrant_valid",
        items = ["rx.radar_valid_out", "rx.radar_out_real_range", "rx.radar_out_imag_range"],
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

    # ──────────────────────────────────────────────────────────────────────
    # report_phase : export data and print summary to terminal
    # ──────────────────────────────────────────────────────────────────────
    def report_phase(self):
        super().report_phase()
        coverage_db.export_to_xml("rx_functional_coverage.xml")
        
        self.logger.info("=========================================================")
        self.logger.info("            RX TOP FUNCTIONAL COVERAGE REPORT            ")
        self.logger.info("=========================================================")
        coverage_db.report_coverage(self.logger.info, bins=True)
        self.logger.info("=========================================================")