"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: monitor.py

    Description:
        This module defines the passive PUVM Monitor for the DDS TX datapath.
        It strictly observes the physical hardware interface without driving 
        any signals, ensuring non-intrusive data collection.

        Key responsibilities:
        1. Passive Sampling: Continuously monitors the DUT pins on active clock edges.
        2. Transaction Reconstruction: Packs the observed pin-level signals back 
           into abstract `seq_item` objects.
        3. Broadcasting: Writes the reconstructed items out through an analysis port 
           so that the Scoreboard and Coverage collectors can independently assess 
           the behavior of the DDS module.
"""

import pyuvm
from pyuvm import *
from seq_item import *
import cocotb
from cocotb.triggers import Timer
from cocotb.clock import Clock
from cocotb.triggers import RisingEdge
class monitor(uvm_monitor):
    def build_phase(self):
        self.mon_port = uvm_analysis_port("mon_port", self)
        self.dut_mon= ConfigDB().get(self,"","DUT")
        
    async def run_phase(self):
        while True:
            await RisingEdge(cocotb.dut_mon.clk)
            
            item = seq_item.create("item")
            # here we will read the signals from the DUT and write them to the analysis port
            
            
            
            self.mon_port.write(item)
            