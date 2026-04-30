"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: driver.py

    Description:
        This module defines the active UVM Driver for the DDS hardware interface.
        It is responsible for translating high-level, abstract data transactions 
        into cycle-accurate, pin-level stimulus on the Device Under Test (DUT).

        Key responsibilities:
        1. Transaction Pulling: Fetches `seq_item` packets from the Sequencer.
        2. Pin Wiggling: Drives specific hardware signals (e.g., FTW, Phase Increment, 
           Enables, and Resets) synchronously with the DDS clock edges.
        3. Protocol Adherence: Ensures all setup, hold, and initialization timings 
           required by the 5G/6G DDS architecture are strictly met.
"""

import cocotb
from cocotb.triggers import RisingEdge, FallingEdge
from pyuvm import *
from dds_seq_item import dds_seq_item

class dds_driver(uvm_driver): 
    def build_phase(self):
        super().build_phase()
        self.dut = ConfigDB().get(self,"","DUT")
          
    async def run_phase(self):
        # 1. Initialize hardware
        self.dut.rst_n.value = 1
        self.dut.enable.value = 0
        self.dut.FTW_start.value = 0
        self.dut.FTW_step.value = 0
        self.dut.cycles.value = 0
        

        for _ in range(5):
            await RisingEdge(self.dut.clk)

        # 2. Main Loop
        while True: 
            req = await self.seq_item_port.get_next_item()
            await FallingEdge(self.dut.clk) 
            
            self.dut.rst_n.value = req.rst_n
            
            if req.rst_n == 0:
                self.dut.enable.value = 0
                await RisingEdge(self.dut.clk)
                self.seq_item_port.item_done()
                continue

            self.dut.FTW_start.value = req.FTW_start
            self.dut.FTW_step.value = req.FTW_step
            self.dut.cycles.value = req.cycles
            self.dut.rst_n.value = req.rst_n
            self.dut.enable.value = 1  
            
            # 3. BLOCK SIMULATION for the chirp duration
            for _ in range(req.cycles):
                await RisingEdge(self.dut.clk)
                
            await FallingEdge(self.dut.clk)
            self.dut.enable.value = 0
            
            self.seq_item_port.item_done()