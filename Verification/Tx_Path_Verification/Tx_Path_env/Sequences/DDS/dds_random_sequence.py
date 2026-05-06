import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_random_seq(uvm_sequence):
    def __init__(self, name="dds_random_seq"):
        super().__init__(name)
        self.num_items = 1  
        
    async def body(self):
        req = dds_seq_item("req")
        req.randomize()
        req.rst_n = 1 # Ensure reset is inactive for random transactions
        req.enable = 1  # Ensure enable is high for valid transactions
        
        
        for _ in range(req.cycles+1):
            
            await self.start_item(req)
            await self.finish_item(req)
       
        req = dds_seq_item("req")
        await self.start_item(req)
        await self.finish_item(req)
        req = dds_seq_item("req")
        await self.start_item(req)
        await self.finish_item(req)