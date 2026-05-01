import pyuvm
from pyuvm import *
from top_seq_item import *  # Importing your updated randomized item
import cocotb
class reset_before_frame_seq(uvm_sequence):
    """
    TC-013: Send frame 1. Insert 10,000 idle cycles. Send frame 2.
    Ensures no stale pipeline state from frame 1 corrupts frame 2.
    """
    def __init__(self, name="reset_seq"):
        super().__init__(name)
        
    async def body(self):
        # ---------------------------------------------------------
        # reset asserted
        # ---------------------------------------------------------
        for _ in range(10):
            req_frame = top_item.create("req_frame")
            await self.start_item(req_frame)
            req_frame.randomize()
            req_frame.rst_n = 0
            await self.finish_item(req_frame)
        
     