import pyuvm
from pyuvm import *
from top_seq_item import top_item
import cocotb
import random
class valid_out_gating_seq(uvm_sequence):
    """
    TC-008: valid_out Gating Verification
    Stimulus: Send exactly one frame, followed by a guaranteed idle period 
    so the monitor can verify valid_out drops strictly to 0 between frames.
    """
    def __init__(self, name="valid_out_gating_seq"):
        super().__init__(name)
        
    async def body(self):
        # 1. Send the exact 4096-cycle frame
        req_frame = top_item.create("req_frame")
        await self.start_item(req_frame)
        
        req_frame.randomize()
        req_frame.rst_n = 1
        req_frame.dds_enable = 1
        req_frame.ofdm_rd_en = 1
        req_frame.cycles = 4096
        
        await self.finish_item(req_frame)
        
        # 2. Force the driver to go idle so the Monitor can verify valid_out goes LOW
        req_idle = top_item.create("req_idle")
        await self.start_item(req_idle)
        
        req_idle.randomize()
        req_idle.rst_n = 1
        req_idle.dds_enable = 0
        req_idle.ofdm_rd_en = 0
        req_idle.cycles = 20  # Wait 20 cycles after the frame finishes
        
        await self.finish_item(req_idle)