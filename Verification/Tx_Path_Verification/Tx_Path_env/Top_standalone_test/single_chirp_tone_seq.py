import pyuvm
from pyuvm import *
from top_seq_item import *  # Importing your updated randomized item
import cocotb

from top_seq_item import top_item # Ensure this is imported to access the set_fixed_chirp function

class single_chirp_tone_seq(uvm_sequence):
    """
    TC-004: Run 100 frames with a single fixed chirp tone (sine wave).
    """
    def __init__(self, name="random_frames_regression_seq"):
        super().__init__(name)
        self.num_frames = 1
        
    async def body(self):
        
        # 2. MAIN FRAME LOOP ('i' is defined here)
        for i in range(self.num_frames):
            
            # Active Data Cycles
            for _ in range(4098):
                req = top_item.create(f"req_frame_{i}")
                await self.start_item(req)
                
                # Force system active
                if _ == 0:
                    req.rst_n      = 1
                    req.dds_enable = 0
                    
                    req.set_fixed_chirp(f0=50e6, B=0e6)  # for sine wave, there is no bandwidth (B=0), f0 is our constant frequency
                else:
                    req.rst_n      = 1
                    req.dds_enable = 1
                    req.set_fixed_chirp(f0=50e6, B=0e6)  # Example fixed chirp for subsequent frames
                    req.cycles     = 4096
                
                await self.finish_item(req)
            
            # Idle Cycles
            for _ in range(4096 * 5):
                req = top_item.create(f"req_idle_{i}")
                await self.start_item(req)
                
                # Stable idle state
                req.rst_n = 1
                req.dds_enable = 0
                req.FTW_start = 0
                req.cycles = 0
                req.FTW_step = 0
                
                await self.finish_item(req)