import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import *
class dds_cyclestress_seq(uvm_sequence):
    def __init__(self, name="dds_cyclestress_seq"):
        super().__init__(name)
    async def body(self):
        # Force the cycle counter to overflow almost instantly
        short_cycles = [1, 2, 3, 5]
        for c in short_cycles:
            req = dds_seq_item("req")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.enable = 1  # Ensure enable is high for valid transactions
            req.cycles = c  # Force short cycle counts to stress the counter logic
            await self.finish_item(req)    