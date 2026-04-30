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
from dds_seq_item import *
import cocotb
from cocotb.clock import Clock
from cocotb.triggers import *

class dds_monitor(uvm_monitor):
    def build_phase(self):
        self.mon_ap = uvm_analysis_port("mon_ap", self)
        self.dut_mon = ConfigDB().get(self,"","DUT")
        
    async def run_phase(self):
        while True:
            # 1. Wait for the clock edge and ensure all signals have settled (ReadOnly)
            await RisingEdge(self.dut_mon.clk)
            await ReadOnly()
            
            # 2. Create a new transaction item to hold the sampled data
            seq_item = dds_seq_item("seq_item")
            
            # 3. Read the input signals from the DUT
            seq_item.rst_n = int(self.dut_mon.rst_n.value)
            seq_item.FTW_start = int(self.dut_mon.FTW_start.value)
            seq_item.FTW_step = int(self.dut_mon.FTW_step.value)
            seq_item.cycles = int(self.dut_mon.cycles.value)
            seq_item.enable = int(self.dut_mon.enable.value)
            
            # 4. Read the output signals from the DUT
            # Note: We use .integer or .signed_integer depending on your SV logic type to cast correctly
            seq_item.final_amplitude = int(self.dut_mon.final_amplitude.value)
            seq_item.valid_out = int(self.dut_mon.valid_out.value)
            
            self.logger.debug(f"Monitoring Item: {seq_item.convert2string_stimulus()}")
            
            # 5. Broadcast the reconstructed transaction to the scoreboard/coverage
            self.mon_ap.write(seq_item)
            
            
            