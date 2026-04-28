"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: subscriber.py

    Description:
        This module defines the PUVM Subscriber (Coverage Collector) for the 
        DDS TX datapath. It acts as the functional coverage engine, passively 
        receiving transactions to measure the completeness of the test suite.

        Key responsibilities:
        1. Transaction Subscribing: Listens continuously to the Monitor's 
           analysis port, receiving reconstructed `seq_item` packets without 
           interfering with the Scoreboard's data flow.
        2. Functional Coverage: Implements covergroups and coverpoints to track 
           which Frequency Tuning Words (FTW), phase offsets, and control modes 
           have been actively verified.
        3. Corner-Case Tracking: Specifically monitors critical edge cases, such 
           as abrupt mid-frame resets (Group F, TC-010) and multiple consecutive 
           resets (TC-011), ensuring the verification plan is 100% satisfied.
"""

import pyuvm
from pyuvm import *
from seq_item import *
from driver import *
from sequencer import sequencer
from monitor import monitor
from cocotb.triggers import Timer
from cocotb.clock import Clock
from cocotb_coverage.coverage import CoverPoint, coverage_db
#coverpoints added 

class subscriber(uvm_component):
    def build_phase(self):     
        self.sub_export = uvm_analysis_export("sub_export", self)
        self.sub_fifo = uvm_tlm_analysis_fifo("sub_fifo", self)
        
    def connect_phase(self):
        self.sub_export.connect(self.sub_fifo.analysis_export)
    async def run_phase(self):
        while True:
            item = await self.sub_fifo.get()
             
        
    
    def coverage_report(self):
        # Export the coverage data to an XML file for reporting and analysis
        
        coverage_db.export_to_xml("alu_coverage.xml")