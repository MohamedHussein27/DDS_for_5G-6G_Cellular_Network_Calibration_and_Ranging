"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_driver.py

    Description:
        The driver translates high-level UVM sequence items into pin-level wiggles.
        It synchronizes with the system clock and drives the valid_in, in_real, 
        and in_imag ports of the fft_4096_top block. It only operates if the 
        parent agent is configured as UVM_ACTIVE.
"""

import cocotb
from cocotb.triggers import *
from pyuvm import *
from ifft_item import *

class ifft_driver(uvm_driver):
    def __init__(self, name, parent):
        super().__init__(name, parent)

    def build_phase(self):
        self.dut = ConfigDB().get(self, "", "DUT")

    async def run_phase(self):

        while True:
            # 1. Get the next transaction from the sequencer
            stim_seq_item = ifft_item("stim_seq_item")
            stim_seq_item = await self.seq_item_port.get_next_item()
            
            # 2. Wait for the clock edge to drive synchronous logic
            await FallingEdge(self.dut.clk)
            
            # 3. Drive the physical pins
            self.dut.rst_n.value    = stim_seq_item.rst_n
            self.dut.valid_in.value = stim_seq_item.valid_in
            self.dut.in_real.value  = stim_seq_item.data_real_in
            self.dut.in_imag.value  = stim_seq_item.data_imag_in
            
            # 4. Notify the sequencer that the transaction is complete
            self.seq_item_port.item_done()