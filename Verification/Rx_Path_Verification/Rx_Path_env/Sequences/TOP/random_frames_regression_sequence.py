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
        self.num_frames = 100
        
    async def body(self):
        for i in range(self.num_frames):
            req = top_item.create(f"req_frame_{i}")
            await self.start_item(req)
            
            # The custom randomize() method inside top_item perfectly 
            # handles the 32-bit registers and 16-bit QAM bounds!
            req.randomize()
            
            # Force system active
            req.rst_n = 1
            req.dds_enable = 1
            req.ofdm_rd_en = 1
            
            # Note: req.cycles is left fully randomized here to test 
            # wildly different frame lengths across the 100 runs!
            
            await self.finish_item(req)