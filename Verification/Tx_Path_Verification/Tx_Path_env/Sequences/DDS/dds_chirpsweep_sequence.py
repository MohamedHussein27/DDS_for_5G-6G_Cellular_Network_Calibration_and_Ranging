import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_chirpsweep_seq(uvm_sequence):
    
    def __init__(self, name="dds_chirpsweep_seq"):
        super().__init__(name)
    async def body(self):
        # Array of aggressive tuning steps
        chirp_steps = [0x100, 0x5000, 0x1FFFF, 0x8FFFF]
        for step_val in chirp_steps:
            req = dds_seq_item("req")
            # Fix the cycles so the wave has time to propagate, force the step
            req.randomize()
            req.rst_n = 1
            req.enable = 1  # Ensure enable is high for valid transactions
            req.cycles = 1000  # Give enough time for the chirp to evolve
            req.FTW_step = step_val  # Force aggressive chirp steps
            for _ in range(req.cycles+1):
                
                await self.start_item(req)
                await self.finish_item(req)
        
            req = dds_seq_item("req")
            await self.start_item(req)
            await self.finish_item(req)
            
        req = dds_seq_item("req")
        await self.start_item(req)
        await self.finish_item(req)
        req = dds_seq_item("req")
        await self.start_item(req)
        await self.finish_item(req)        