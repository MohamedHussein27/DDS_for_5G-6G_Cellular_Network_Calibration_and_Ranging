import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_fft_boundary_seq(uvm_sequence):
    def __init__(self, name="dds_fft_boundary_seq"):
        super().__init__(name)

    async def body(self):
        # 4095 = 1 cycle short (Underflow)
        # 4096 = Perfect FFT Symbol
        # 4097 = 1 cycle over (Overflow/Backpressure)
        boundaries = [4095, 4096, 4097]
        
        for c in boundaries:
            req = dds_seq_item("req")
            
            req.randomize()
            req.enable = 1  # Ensure enable is high for valid transactions
            req.rst_n = 1
            req.cycles = c  # Force the specific FFT boundary
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