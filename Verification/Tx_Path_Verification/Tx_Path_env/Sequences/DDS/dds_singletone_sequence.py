import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_singletone_seq(uvm_sequence):
    def __init__(self, name="dds_singletone_seq"):
        super().__init__(name)
        self.num_items = 15
    async def body(self):
        
        # Generate 15 pure single tones (No chirp acceleration)
        for _ in range(self.num_items):
            req = dds_seq_item("req")
            await self.start_item(req)
            # Force FTW_step to 0, randomize the rest
            req.randomize()
            req.rst_n = 1
            req.FTW_step = 0
            req.enable = 1  # Ensure enable is high for valid transactions
            await self.finish_item(req)
