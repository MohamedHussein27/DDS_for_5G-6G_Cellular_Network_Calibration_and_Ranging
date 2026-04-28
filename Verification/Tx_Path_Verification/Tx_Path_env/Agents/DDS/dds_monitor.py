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
        self.mon_port = uvm_analysis_port("mon_port", self)
        self.dut_mon = ConfigDB().get(self,"","DUT")
        
    async def run_phase(self):
        while True:
            # 1. Wait for the clock edge and ensure all signals have settled (ReadOnly)
            await RisingEdge(self.dut_mon.clk)
            await ReadOnly()
            
            # 2. Create a new transaction item to hold the sampled data
            seq_item = dds_seq_item("seq_item")
            
            # 3. Read the input signals from the DUT
            seq_item.rst_n = self.dut_mon.rst_n.value.integer
            seq_item.FTW_start = self.dut_mon.FTW_start.value.integer
            seq_item.FTW_step = self.dut_mon.FTW_step.value.integer
            seq_item.cycles = self.dut_mon.cycles.value.integer
            
            # 4. Read the output signals from the DUT
            # Note: We use .integer or .signed_integer depending on your SV logic type to cast correctly
            seq_item.final_amplitude = self.dut_mon.final_amplitude.value.integer 
            
            self.logger.debug(f"Monitoring Item: {seq_item.convert2string_stimulus()}")
            
            # 5. Broadcast the reconstructed transaction to the scoreboard/coverage
            self.mon_port.write(seq_item)