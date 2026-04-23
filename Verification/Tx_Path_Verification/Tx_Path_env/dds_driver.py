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

import pyuvm
from pyuvm import *
from seq_item import *
import cocotb
from cocotb.triggers import Timer
from cocotb.clock import Clock
from cocotb.triggers import FallingEdge

class driver(uvm_driver): 
    
    def build_phase(self):
        self.dut_drv = ConfigDB().get(self,"","DUT")
          
    async def run_phase(self):
        while True: 
            seq_item = seq_item.create("seq_item")
            seq_item = await self.seq_item_port.get_next_item()
            await FallingEdge(self.dut_drv.clk) 
            
            # here we will drive the signals to the DUT using the values from the sequence item
            
            
            
            self.seq_item_port.item_done()