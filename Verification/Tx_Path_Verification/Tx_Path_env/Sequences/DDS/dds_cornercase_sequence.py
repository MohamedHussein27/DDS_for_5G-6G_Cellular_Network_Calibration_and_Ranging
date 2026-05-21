import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_cornercase_seq(uvm_sequence):
    def __init__(self, name="dds_cornercase_seq"):
        super().__init__(name)

    async def body(self):
        # Test Absolute Maximums 
        req = dds_seq_item("req")
        await self.start_item(req)
        req.randomize()
        req.enable = 1  # Ensure enable is high for valid transactions
        req.rst_n = 1
        req.FTW_start = 0xFFFFFFFF  # All 1s
        req.FTW_step = 0xFFFFFFFF   # All 1s
        req.cycles = 0x1FFF         # Max 13-bit value (8191)
        
        await self.finish_item(req)