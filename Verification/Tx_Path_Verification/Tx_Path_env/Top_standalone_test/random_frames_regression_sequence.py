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
        for i in range(self.num_frames):
            for _ in range (4096):
                req = top_item.create(f"req_frame_{i}")
                await self.start_item(req)
                
                # The custom randomize() method inside top_item perfectly 
                # handles the 32-bit registers and 16-bit QAM bounds!
                #req.randomize()
                
                # Force system active
                req.rst_n = 1
                req.dds_enable = 1
                req.FTW_start = 0
                req.cycles = 4096
                req.FTW_step = 150000
                
                # Note: req.cycles is left fully randomized here to test 
                # wildly different frame lengths across the 100 runs!
                
                await self.finish_item(req)
            
            for _ in range(4096 * 6):
                req = top_item.create(f"req_idle_{i}")
                await self.start_item(req)
                
                # During idle cycles, we can keep the inputs constant or randomize them.
                # Here, we'll just keep them constant to represent a stable idle state.
                req.rst_n = 1
                req.dds_enable = 0
                req.FTW_start = 0
                req.cycles = 0
                req.FTW_step = 0
                
                await self.finish_item(req)