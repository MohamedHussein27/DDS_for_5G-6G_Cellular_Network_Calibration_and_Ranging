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
 
from top_seq_item import top_item
 
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
    # write : The mandatory uvm_subscriber hook (Replaces run_phase/FIFO)
    # ──────────────────────────────────────────────────────────────────────
    def write(self, item):
        """
        This function is instantly called by the Monitor's Analysis Port. 
        It executes in zero simulation time.
        """
        # 1. Skip backdoor transactions (they are purely software constructs)
        if hasattr(item, 'is_backdoor') and item.is_backdoor:
            return
            
        # 2. Pass the transaction to the coverage sampling decorators
        self.sample(item)

 
    # ══════════════════════════════════════════════════════════════════════
    # COVER POINTS  (all bins defined with lambda pin functions)
    # ══════════════════════════════════════════════════════════════════════

    # ── CP1 : wr_en ───────────────────────────────────────────────────────
    @CoverPoint(
        "top.wr_en",
        xf   = lambda tr: int(tr.wr_en),
        bins = [0, 1],
        bins_labels = ["wr_inactive", "wr_active"],
    )
    # ── CP2 : rd_en ───────────────────────────────────────────────────────
    @CoverPoint(
        "top.rd_en",
        xf   = lambda tr: int(tr.rd_en),
        bins = [0, 1],
        bins_labels = ["rd_inactive", "rd_active"],
    )
    # ── CP3 : addr ────────────────────────────────────────────────────────
    @CoverPoint(
        "top.addr",
        xf   = lambda tr: int(tr.addr),
        bins = [0x0, 0x4, 0x8],
        bins_labels = ["ADDR_START", "ADDR_STEP", "ADDR_CYCLES"],
    )
    # ── CP4 : dds_ready_flag ──────────────────────────────────────────────
    @CoverPoint(
        "top.dds_ready_flag",
        xf   = lambda tr: int(tr.dds_ready_flag),
        bins = [0, 1],
        bins_labels = ["dds_not_ready", "dds_ready"],
    )
    # ── CP5 : tx_valid ────────────────────────────────────────────────────
    @CoverPoint(
        "top.tx_valid",
        xf   = lambda tr: int(tr.tx_valid),
        bins = [0, 1],
        bins_labels = ["tx_valid_inactive", "tx_valid_active"],
    )
    # ── CP6 : rst_n ───────────────────────────────────────────────────────
    @CoverPoint(
        "top.rst_n",
        xf   = lambda tr: int(tr.rst_n),
        bins = [0, 1],
        bins_labels = ["reset_active", "reset_inactive"],
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
    
    # ══════════════════════════════════════════════════════════════════════
    # CROSS COVERAGE
    # ══════════════════════════════════════════════════════════════════════
    
    # ── CP9 : Ensure we WROTE to every address ────────────────────────────
    @CoverCross(
        "top.bus_write_cross",
        items = ["top.addr", "top.wr_en"],
        # We only care about ensuring wr_en=1 happened for 0x0, 0x4, 0x8
        ign_bins = [("ADDR_START", "wr_inactive"), 
                    ("ADDR_STEP", "wr_inactive"), 
                    ("ADDR_CYCLES", "wr_inactive")]
    )
    
    # ── CP10 : Ensure we READ from every address ──────────────────────────
    @CoverCross(
        "top.bus_read_cross",
        items = ["top.addr", "top.rd_en"],
        # We only care about ensuring rd_en=1 happened for 0x0, 0x4, 0x8
        ign_bins = [("ADDR_START", "rd_inactive"), 
                    ("ADDR_STEP", "rd_inactive"), 
                    ("ADDR_CYCLES", "rd_inactive")]
    )

    # ── CP11 : tx_out_real × tx_out_imag quadrant, tx_valid cross ─────────
    @CoverCross(
        "top.tx_output_quadrant_valid",
        items = ["top.tx_valid", "top.tx_out_real_range", "top.tx_out_imag_range"],
    )

    def sample(self, tr):
        pass # Used by decorators
 
    # ──────────────────────────────────────────────────────────────────────
    # report_phase : Generate the XML and ASCII coverage summary
    # ──────────────────────────────────────────────────────────────────────
    def report_phase(self):
        # 1. Export the detailed coverage to an XML file for deeper inspection later
        xml_file = "tx_top_coverage.xml"
        coverage_db.export_to_xml(xml_file)

        # 2. Build the custom ASCII terminal report
        self.logger.info("╔══════════════════════════════════════════════════════════════╗")
        self.logger.info("║               TX_TOP FUNCTIONAL COVERAGE REPORT              ║")
        self.logger.info("╠══════════════════════════════════════════════════════════════╣")

        all_pass = True

        # Loop through every coverpoint/cross defined in the coverage database
        for name in coverage_db:
            # Only print coverage items belonging to the top-level block
            # (Assuming your coverpoints were defined with "top." like "top.wr_en")
            if name.startswith("top.") and not name.startswith("top.ifft"):
                cp = coverage_db[name]
                
                # Retrieve the coverage percentage (0.0 to 100.0)
                percentage = cp.cover_percentage
                
                # Check if it reached the 100% goal
                if percentage >= 100.0:
                    status = "PASS ✓"
                else:
                    status = "FAIL ✗"
                    all_pass = False

                # Format the line to keep columns perfectly aligned
                # Extracts just the coverpoint name (e.g., 'wr_en' instead of 'top.wr_en')
                short_name = name.split(".")[-1]
                
                self.logger.info(f"║  {short_name:<34} : {status:<8} ({percentage:>5.1f}%)  ║")

        # Overall Status
        overall_status = "PASS ✓" if all_pass else "FAIL ✗"
        self.logger.info("║ ──────────────────────────────────────────────────────────── ║")
        self.logger.info(f"║  OVERALL COVERAGE TARGETS           : {overall_status:<17}  ║")
        self.logger.info("╚══════════════════════════════════════════════════════════════╝")
        
        self.logger.info(f"Detailed coverage data exported to {xml_file}")