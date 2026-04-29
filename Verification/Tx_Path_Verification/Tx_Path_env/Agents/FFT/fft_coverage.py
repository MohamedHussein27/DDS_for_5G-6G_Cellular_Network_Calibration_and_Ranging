"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: fft_subscriber.py

    Description:
        This functional coverage subscriber passively monitors the FFT boundary 
        to ensure the verification environment adequately stimulates the datapath. 
        It utilizes cocotb_coverage to track control signals (valid_in, valid_out, 
        rst_n) and ensures the 16-bit signed fixed-point data extensively hits 
        critical boundary conditions (min_neg, zero, max_pos). It also cross-
        references the real and imaginary components to ensure all four complex 
        quadrants are actively stimulated during valid transmission cycles.
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

from fft_seq_item import fft_item

# ─────────────────────────────────────────────────────────────────────────────
# Constants that mirror the RTL parameters
# ─────────────────────────────────────────────────────────────────────────────
WL   = 16                          # word length (signed)
N    = 4096                        # FFT size
MAX_VAL  =  (1 << (WL - 1)) - 1    #  32767
MIN_VAL  = -(1 << (WL - 1))        # -32768
MID_POS  =  MAX_VAL // 2           #  16383  (upper half of positive range)
MID_NEG  =  MIN_VAL // 2           # -16384  (lower half of negative range)

# ─────────────────────────────────────────────────────────────────────────────
# Helper – sign-extend a raw integer read from a 16-bit RTL signal
# ─────────────────────────────────────────────────────────────────────────────
def _sign16(v):
    v = int(v) & 0xFFFF
    if v & 0x8000:
        return v - 0x10000
    else:
        return v


class fft_subscriber(uvm_subscriber):
    """
    Functional-coverage subscriber for the 4096-point FFT.

    Cover points
    ─────────────
    CP1  valid_in        : both de-asserted (0) and asserted (1)
    CP2  valid_out       : both de-asserted (0) and asserted (1)
    CP3  in_real range   : negative / zero / positive / max-positive / min-negative
    CP4  in_imag range   : negative / zero / positive / max-positive / min-negative
    CP5  out_real range  : negative / zero / positive / max-positive / min-negative
    CP6  out_imag range  : negative / zero / positive / max-positive / min-negative
    CP7  rst_n           : reset active (0) and inactive (1)
    CP8  in_real × in_imag quadrant cross : (pos,pos) (pos,neg) (neg,pos) (neg,neg)
    CP9  valid_out × quadrant cross       : Ensures valid data crossed all quadrants

    All bins are implemented with lambda "pin" functions so the coverage
    decorator never needs to reference the class directly.
    """

    def __init__(self, name, parent):
        super().__init__(name, parent)

    # ──────────────────────────────────────────────────────────────────────
    # build_phase : create export + FIFO (mirroring the sample template)
    # ──────────────────────────────────────────────────────────────────────
    def build_phase(self):
        self.sub_export = uvm_analysis_export("sub_export", self)
        self.sub_fifo   = uvm_tlm_analysis_fifo("sub_fifo", self)
        # re-point the export straight at the FIFO's analysis_export
        self.sub_export = self.sub_fifo.analysis_export

    # ──────────────────────────────────────────────────────────────────────
    # write : mandatory uvm_subscriber hook
    # ──────────────────────────────────────────────────────────────────────
    def write(self, item):
        self.logger.debug(f"FFTSubscriber received: {item.convert2string()}")

    # ══════════════════════════════════════════════════════════════════════
    # COVER POINTS  (all bins defined with lambda pin functions)
    # ══════════════════════════════════════════════════════════════════════

    # ── CP1 : valid_in ────────────────────────────────────────────────────
    @CoverPoint(
        "top.fft.valid_in",
        xf   = lambda tr: int(tr.valid_in),
        bins = [0, 1],
        bins_labels = ["valid_in_deasserted", "valid_in_asserted"],
    )
    # ── CP2 : valid_out ───────────────────────────────────────────────────
    @CoverPoint(
        "top.fft.valid_out",
        xf   = lambda tr: int(tr.valid_out),
        bins = [0, 1],
        bins_labels = ["valid_out_deasserted", "valid_out_asserted"],
    )
    # ── CP3 : in_real range ───────────────────────────────────────────────
    @CoverPoint(
        "top.fft.in_real_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.in_real) == MIN_VAL  else
            "negative" if _sign16(tr.in_real) <  0        else
            "zero"     if _sign16(tr.in_real) == 0        else
            "max_pos"  if _sign16(tr.in_real) == MAX_VAL  else
            "positive"
        ),
        bins        = ["min_neg", "negative", "zero", "positive", "max_pos"],
        bins_labels = ["in_real_min_neg", "in_real_neg",
                       "in_real_zero", "in_real_pos", "in_real_max_pos"],
    )
    # ── CP4 : in_imag range ───────────────────────────────────────────────
    @CoverPoint(
        "top.fft.in_imag_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.in_imag) == MIN_VAL  else
            "negative" if _sign16(tr.in_imag) <  0        else
            "zero"     if _sign16(tr.in_imag) == 0        else
            "max_pos"  if _sign16(tr.in_imag) == MAX_VAL  else
            "positive"
        ),
        bins        = ["min_neg", "negative", "zero", "positive", "max_pos"],
        bins_labels = ["in_imag_min_neg", "in_imag_neg",
                       "in_imag_zero", "in_imag_pos", "in_imag_max_pos"],
    )
    # ── CP5 : out_real range ──────────────────────────────────────────────
    @CoverPoint(
        "top.fft.out_real_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.out_real) == MIN_VAL  else
            "negative" if _sign16(tr.out_real) <  0        else
            "zero"     if _sign16(tr.out_real) == 0        else
            "max_pos"  if _sign16(tr.out_real) == MAX_VAL  else
            "positive"
        ),
        bins        = ["min_neg", "negative", "zero", "positive", "max_pos"],
        bins_labels = ["out_real_min_neg", "out_real_neg",
                       "out_real_zero", "out_real_pos", "out_real_max_pos"],
    )
    # ── CP6 : out_imag range ──────────────────────────────────────────────
    @CoverPoint(
        "top.fft.out_imag_range",
        xf = lambda tr: (
            "min_neg"  if _sign16(tr.out_imag) == MIN_VAL  else
            "negative" if _sign16(tr.out_imag) <  0        else
            "zero"     if _sign16(tr.out_imag) == 0        else
            "max_pos"  if _sign16(tr.out_imag) == MAX_VAL  else
            "positive"
        ),
        bins        = ["min_neg", "negative", "zero", "positive", "max_pos"],
        bins_labels = ["out_imag_min_neg", "out_imag_neg",
                       "out_imag_zero", "out_imag_pos", "out_imag_max_pos"],
    )
    # ── CP7 : rst_n ───────────────────────────────────────────────────────
    @CoverPoint(
        "top.fft.rst_n",
        xf   = lambda tr: int(tr.rst_n),
        bins = [0, 1],
        bins_labels = ["reset_active", "reset_inactive"],
    )
    # ── CP8 : in_real × in_imag quadrant cross ───────────────────────────
    @CoverCross(
        "top.fft.input_quadrant",
        items = ["top.fft.in_real_range", "top.fft.in_imag_range"],
    )
    # ── CP9 : in_real × in_imag quadrant, valid out cross ────────────────
    @CoverCross(
        "top.fft.input_quadrant_valid_out",
        items = ["top.fft.valid_out", "top.fft.in_real_range", "top.fft.in_imag_range"],
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