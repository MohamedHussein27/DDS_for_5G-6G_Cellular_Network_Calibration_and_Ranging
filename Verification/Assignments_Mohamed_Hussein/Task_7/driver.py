import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

from sequence_item import *

class Driver(uvm_driver):

    def __init__(self, name, parent):
        super().__init__(name, parent)

    def build_phase(self):
        self.dut_drv = ConfigDB().get(self,"","DUT")
    

    async def run_phase(self):
        while True:
            stim_seq_item = Sequence_Item.create("stim_seq_item")
            stim_seq_item = await self.seq_item_port.get_next_item()

            await FallingEdge(self.dut_drv.clk)
            # driving
            # Displaying the values
            """self.logger.info(f"[DRIVER] sending: rst_n={stim_seq_item.rst_n}, "
                         f"a={stim_seq_item.a}, b={stim_seq_item.b}, op={stim_seq_item.op}")"""

            self.dut_drv.rst_n.value = stim_seq_item.rst_n
            self.dut_drv.a.value     = stim_seq_item.a
            self.dut_drv.b.value     = stim_seq_item.b
            self.dut_drv.op.value    = stim_seq_item.op

            self.seq_item_port.item_done()
  