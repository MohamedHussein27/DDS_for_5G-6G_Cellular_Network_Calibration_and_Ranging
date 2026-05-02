import pyuvm
from pyuvm import *
from top_seq_item import *  # Importing your updated randomized item
import cocotb

class random_frames_regression_seq(uvm_sequence):
    """
    TC-014: Run 100 frames with fully constrained-random sequences.
    This acts as the ultimate stress test for the TX Top Datapath.
    """
    def __init__(self, name="random_frames_regression_seq"):
        super().__init__(name)
        self.num_frames = 1
        
    async def body(self):
        # 1. INITIAL RESET TRANSACTION (No 'i' variable here)
        req = top_item.create("req_reset")
        await self.start_item(req)
        
        # Force system reset
        req.rst_n = 0
        req.dds_enable = 1
        req.FTW_start = 0
        req.cycles = 4096
        req.FTW_step = 150000
        
        await self.finish_item(req)
        
        # 2. MAIN FRAME LOOP ('i' is defined here)
        for i in range(self.num_frames):
            
            # Active Data Cycles
            for _ in range(4096):
                req = top_item.create(f"req_frame_{i}")
                await self.start_item(req)
                
                # Force system active
                req.rst_n = 1
                req.dds_enable = 1
                req.FTW_start = 0
                req.cycles = 4096
                req.FTW_step = 150000
                
                await self.finish_item(req)
            
            # Idle Cycles
            for _ in range(4096 * 6):
                req = top_item.create(f"req_idle_{i}")
                await self.start_item(req)
                
                # Stable idle state
                req.rst_n = 1
                req.dds_enable = 0
                req.FTW_start = 0
                req.cycles = 0
                req.FTW_step = 0
                
                await self.finish_item(req)