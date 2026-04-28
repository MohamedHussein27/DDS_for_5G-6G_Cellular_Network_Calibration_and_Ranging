"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: scoreboard.py

    Description:
        This module defines the PUVM Scoreboard for the DDS TX datapath.
        It serves as the ultimate verification authority, determining the 
        pass/fail status of the DUT by comparing actual hardware outputs 
        against an ideal, bit-accurate reference model.

        Key responsibilities:
        1. Reference Modeling: Calculates expected sine/cosine waveforms mathematically 
           based on the applied Frequency Tuning Words (FTW), phase increments, 
           and control signals.
        2. Data Collection: Captures observed `seq_item` transactions from the 
           passive Monitors via UVM analysis exports.
        3. Comparison & Tolerance: Compares the hardware's generated digital 
           waveform against the golden model, accounting for expected 
           quantization noise and phase truncation inherent to DDS architecture.
        4. Reporting: Logs precise mismatch locations, tracks verification 
           statistics, and flags fatal errors if functional requirements are violated.
"""
import pyuvm
from pyuvm import *
from dds_seq_item import *
from dds_driver import dds_driver
from dds_monitor import dds_monitor
from dds_sequencer import dds_sequencer
from cocotb.triggers import Timer
from cocotb.clock import Clock

class dds_scoreboard(uvm_scoreboard):
    def build_phase(self):
        # 1. Initialize counters for the += operation to work
        self.passed_test_cases = 0
        self.failed_test_cases = 0
        
        self.sc_export = uvm_analysis_export("sc_export", self)
        self.sc_fifo = uvm_tlm_analysis_fifo("sc_fifo", self)
        
    def connect_phase(self):
        self.sc_export.connect(self.sc_fifo.analysis_export)
        
    # 2. Add async here
    async def run_phase(self):
        while True:
            # 3. Add await here
            item = await self.sc_fifo.get()
            """************************** TEST CASES **************************"""
    def report_phase(self):
        cocotb.log.info("========================================")
        cocotb.log.info("         SCOREBOARD FINAL REPORT        ")
        cocotb.log.info("========================================")
        cocotb.log.info(f" PASSED: {self.passed_test_cases}")
        cocotb.log.info(f" FAILED: {self.failed_test_cases}")
        cocotb.log.info("========================================")
        if self.failed_test_cases == 0:
            cocotb.log.info(" TEST STATUS: SUCCESS ")
        else:
            cocotb.log.error(" TEST STATUS: FAILED ")        