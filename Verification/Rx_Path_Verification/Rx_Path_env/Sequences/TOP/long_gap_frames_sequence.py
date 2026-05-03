import pyuvm
from pyuvm import *
from top_seq_item import *  # Importing your updated randomized item
import cocotb
class long_gap_frames_seq(uvm_sequence):
    """
    TC-013: Send frame 1. Insert 10,000 idle cycles. Send frame 2.
    Ensures no stale pipeline state from frame 1 corrupts frame 2.
    """
    def __init__(self, name="long_gap_frames_seq"):
        super().__init__(name)
        
    async def body(self):
        # ---------------------------------------------------------
        # FRAME 1: Valid Data
        # ---------------------------------------------------------
        req_frame1 = top_item.create("req_frame1")
        await self.start_item(req_frame1)
        req_frame1.randomize()
        req_frame1.rst_n = 1
        req_frame1.dds_enable = 1
        req_frame1.ofdm_rd_en = 1
        req_frame1.cycles = 4096
        await self.finish_item(req_frame1)
        
        # ---------------------------------------------------------
        # THE GAP: 10,000 Idle Cycles
        # ---------------------------------------------------------
        req_idle = top_item.create("req_idle")
        await self.start_item(req_idle)
        req_idle.randomize()
        req_idle.rst_n = 1
        req_idle.dds_enable = 0  # Disable DDS
        req_idle.ofdm_rd_en = 0  # Disable OFDM read
        req_idle.cycles = 10000  # Force driver to wait 10,000 cycles
        await self.finish_item(req_idle)

        # ---------------------------------------------------------
        # FRAME 2: Valid Data (Different Stimulus)
        # ---------------------------------------------------------
        req_frame2 = top_item.create("req_frame2")
        await self.start_item(req_frame2)
        req_frame2.randomize()
        req_frame2.rst_n = 1
        req_frame2.dds_enable = 1
        req_frame2.ofdm_rd_en = 1
        req_frame2.cycles = 4096
        await self.finish_item(req_frame2)