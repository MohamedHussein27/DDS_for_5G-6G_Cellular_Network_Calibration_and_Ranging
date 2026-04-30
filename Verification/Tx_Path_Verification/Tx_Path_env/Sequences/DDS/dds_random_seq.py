import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_random_seq(uvm_sequence):
    def __init__(self, name="dds_random_seq"):
        super().__init__(name)
        # Default to 10, but allow the Test to overwrite this!
        self.num_items = 10 
        
    async def body(self):
        for _ in range(self.num_items):
            req = dds_seq_item("req")
            await self.start_item(req)
            req.randomize()
            req.enable = 1  # Ensure enable is high for valid transactions
            await self.finish_item(req)