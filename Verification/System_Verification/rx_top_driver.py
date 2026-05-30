"""
    Module: rx_driver.py

    Description:
        The RX driver receives rx_item transactions from the sequencer and 
        drives the physical pins of the RX_TOP module. It handles both the 
        channel inputs and the reference RAM inputs simultaneously.
"""

import pyuvm
from pyuvm import *
import cocotb
from cocotb.triggers import FallingEdge

# Assuming you named your sequence item class rx_item
from rx_top_seq_item import *

class rx_driver(uvm_driver):
    def build_phase(self):
        # Fetch the DUT from the UVM Configuration Database
        self.dut_drv = ConfigDB().get(self, "", "DUT")
        
    async def run_phase(self):
        while True:
            # 1. Get the next transaction item from the sequencer
            seq_item = await self.seq_item_port.get_next_item()
            
            # 2. Synchronize with the clock edge before driving the signals.
            # Driving on the falling edge ensures setup and hold times are 
            # perfectly met for the RTL's posedge logic.
            await FallingEdge(self.dut_drv.clk)
            
            # 3. Drive the stimulus into the DUT inputs
            # Global
            
            self.dut_drv.rst_n.value = seq_item.rst_n
            
            # Channel Input Path 
            self.dut_drv.rx_valid_in.value = seq_item.rx_valid_in
            self.dut_drv.rx_in_re.value  = seq_item.rx_in_re
            self.dut_drv.rx_in_im.value  = seq_item.rx_in_im
            
            # Reference RAM Input Path (From TX)
            self.dut_drv.ref_wr_en.value = seq_item.ref_wr_en
            self.dut_drv.ref_wr_re.value = seq_item.ref_wr_re
            self.dut_drv.ref_wr_im.value = seq_item.ref_wr_im
            
            # 4. Notify the sequencer that this transaction is complete
            self.seq_item_port.item_done()