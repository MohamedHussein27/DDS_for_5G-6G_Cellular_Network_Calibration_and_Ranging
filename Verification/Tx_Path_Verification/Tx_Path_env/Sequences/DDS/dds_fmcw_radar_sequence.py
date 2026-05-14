import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_fmcw_radar_seq(uvm_sequence):
    """Generates continuous back-to-back linear chirps for radar sensing."""
    def __init__(self, name="dds_fmcw_radar_seq"):
        super().__init__(name)

    async def body(self):
        for _ in range(4): # Send 4 back-to-back radar pulses
            # 1. The Ramp (Chirp)
            req = dds_seq_item.create("req")
            
            req.rst_n = 1
            req.enable = 1
            req.FTW_start = 0x0
            req.FTW_step = 0x00050000 # Moderate acceleration
            req.cycles = 1024         # Standard FMCW symbol length
            for _ in range(req.cycles+1):
                await self.start_item(req)
                await self.finish_item(req)
            
            # 2. The Flyback (Reset the wave)
            req_flyback = dds_seq_item.create("req_flyback")
            await self.start_item(req_flyback)
            req_flyback.rst_n = 1
            req_flyback.enable = 0    # Pull low to flush pipeline
            req_flyback.cycles = 2    # Brief 2-cycle wait
            await self.finish_item(req_flyback)  
        req_stop = dds_seq_item.create("req_stop")
        await self.start_item(req_stop)
        req_stop.rst_n = 1
        req_stop.enable = 0
        req_stop.cycles = 1
        await self.finish_item(req_stop)      
        req_stop = dds_seq_item.create("req_stop")
        await self.start_item(req_stop)
        req_stop.rst_n = 1
        req_stop.enable = 0
        req_stop.cycles = 1
        await self.finish_item(req_stop)     