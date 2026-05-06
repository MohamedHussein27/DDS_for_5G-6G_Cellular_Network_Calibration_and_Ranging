import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_max_cov_seq(uvm_sequence):
    """A directed sequence designed specifically to hit the FTW_start max bin."""
    def __init__(self, name="dds_max_cov_seq"):
        super().__init__(name)

    async def body(self):
        # 1. Target the absolute "max" bin for FTW_start
        req = dds_seq_item.create("req")
        req.randomize()
        req.rst_n = 1
        req.enable = 1
        req.FTW_start = 0xFFFFFFFF  # <--- This is the golden ticket!
        req.FTW_step = 0x0          # Keep it a simple single tone
        req.cycles = 50             # Just a short 50-cycle burst
        for _ in range(req.cycles+1):
            await self.start_item(req)
            await self.finish_item(req)
            
            
        req2 = dds_seq_item.create("req2")
        
        req2.randomize()
        req2.rst_n = 1
        req2.enable = 1
        req2.FTW_start = 0xFFFFFFFF  # Max Phase
        req2.FTW_step = 0x00050000   # Step > 0 (Chirping!)
        req2.cycles = 50             
        for _ in range(req2.cycles+1):
            await self.start_item(req2)
            await self.finish_item(req2)   
       
        
        # 2. Safely flush and turn off the pipeline
        req_stop = dds_seq_item.create("req_stop")
        await self.start_item(req_stop)
        req_stop.rst_n = 1
        req_stop.enable = 0 
        req_stop.cycles = 2         # Brief wait to register the turn-off
        await self.finish_item(req_stop)
        