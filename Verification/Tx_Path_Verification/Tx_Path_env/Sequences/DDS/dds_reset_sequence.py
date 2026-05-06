import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 

class dds_reset(uvm_sequence):
    def __init__(self, name="dds_reset"):
        super().__init__(name) 
    async def body(self):
            req = dds_seq_item("req")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 0  # FORCE RESET ACTIVE
            req.enable=1
            await self.finish_item(req)    
            req = dds_seq_item("req")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 0  # FORCE RESET ACTIVE
            req.enable=0
            await self.finish_item(req)  
            await self.start_item(req)
            await self.finish_item(req)