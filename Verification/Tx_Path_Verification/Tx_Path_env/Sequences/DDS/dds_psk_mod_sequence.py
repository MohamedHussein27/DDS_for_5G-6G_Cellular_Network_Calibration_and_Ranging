import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_psk_modulation_seq(uvm_sequence):
    """Simulates QPSK symbol transmission by shifting the phase 90 degrees."""
    def __init__(self, name="dds_psk_modulation_seq"):
        super().__init__(name)

    async def body(self):
        # 0°, 90°, 180°, 270° in a 32-bit accumulator
        qpsk_phases = [0x00000000, 0x40000000, 0x80000000, 0xC0000000]
        
        for phase in qpsk_phases:
            req = dds_seq_item.create("req")
            req.randomize()
            req.rst_n = 1
            req.enable = 1
            req.FTW_start = phase
            req.FTW_step = 0x00010000 # Constant carrier frequency
            req.cycles = 50           # Short 50-cycle symbol duration
            for _ in range(req.cycles+1):
                await self.start_item(req)
                await self.finish_item(req)
            req = dds_seq_item.create("req")
            await self.start_item(req)
            await self.finish_item(req)
        req = dds_seq_item("req")
        await self.start_item(req)
        await self.finish_item(req)
        req = dds_seq_item("req")
        await self.start_item(req)
        await self.finish_item(req) 