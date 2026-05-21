"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: subscriber.py

    Description:
        This module defines the PUVM Subscriber (Coverage Collector) for the 
        DDS TX datapath. It acts as the functional coverage engine, passively 
        receiving transactions to measure the completeness of the test suite.
"""

import pyuvm
from pyuvm import *
import cocotb
from cocotb_coverage.coverage import CoverPoint, CoverCross, coverage_db
from dds_seq_item import dds_seq_item  

# =======================================================================
# Coverage Definition (Outside the class, using decorators)
# =======================================================================

# CoverPoint for the Reset Signal
@CoverPoint(
    "dds_coverage.rst_n", 
    xf = lambda item: item.rst_n, 
    bins = [0, 1],
    bins_labels = ["reset_active", "reset_inactive"]
)
# CoverPoint for the Enable Signal
@CoverPoint(
    "dds_coverage.enable", 
    xf = lambda item: item.enable, 
    bins = [0, 1],
    bins_labels = ["disabled", "enabled"]
)
# CoverPoint for FTW_start (Categorized)
@CoverPoint(
    "dds_coverage.FTW_start", 
    xf = lambda item: (
        "zero" if item.FTW_start == 0 else 
        "max" if item.FTW_start == 0xFFFFFFFF else 
        "mid"
    ), 
    bins = ["zero", "mid", "max"]
)
# CoverPoint for FTW_step (Categorized to detect single tone vs chirp)
@CoverPoint(
    "dds_coverage.FTW_step", 
    xf = lambda item: "single_tone" if item.FTW_step == 0 else "chirping", 
    bins = ["single_tone", "chirping"]
)

#  Chirp Direction (Up-chirp vs Down-chirp) 
@CoverPoint(
    "dds_coverage.chirp_direction", 
    xf = lambda item: (
        "single_tone" if item.FTW_step == 0 else 
        "down_chirp"  if item.FTW_step >= 0x80000000 else 
        "up_chirp"
    ), 
    bins = ["single_tone", "up_chirp", "down_chirp"]
)

#  Burst Duration / Cycle Count 
@CoverPoint(
    "dds_coverage.cycles_duration",
    xf = lambda item: (
        "short_burst" if item.cycles <= 10 else
        "fft_boundary" if item.cycles in [4095, 4096, 4097] else
        "long_stream"
    ),
    bins = ["short_burst", "fft_boundary", "long_stream"]
)

#  CROSS COVERAGE: Tone Type crossed with Starting Phase 
@CoverCross(
    "dds_coverage.start_vs_step",
    items = ["dds_coverage.FTW_start", "dds_coverage.FTW_step"]
)

def sample_dds_coverage(item):
    """ This empty function serves as the trigger point for the coverage decorators. """
    pass


# =======================================================================
# Subscriber Component
# =======================================================================
class dds_subscriber(uvm_component):  
    def build_phase(self):     
        super().build_phase()
        
        # 1. Create the FIFO
        self.sub_fifo = uvm_tlm_analysis_fifo("sub_fifo", self)
        
        # 2. Connect the analysis export to the FIFO
        self.analysis_export = self.sub_fifo.analysis_export
    async def run_phase(self):
        while True:
            # Grab the latest transaction from the Monitor's FIFO
            item = await self.sub_fifo.get()
            
            # Trigger the coverage collection
            sample_dds_coverage(item)

    def report_phase(self):
     # 1. Export the XML for external tools
                coverage_file = "dds_coverage.xml"
                coverage_db.export_to_xml(coverage_file)
 
    # 2. Print beautiful coverage tables to the terminal!
                cocotb.log.info("=================================================")
                cocotb.log.info(" DDS COVERAGE REPORT ")
                cocotb.log.info("=================================================")
                coverage_db.report_coverage(cocotb.log.info, bins=True)
                   
               