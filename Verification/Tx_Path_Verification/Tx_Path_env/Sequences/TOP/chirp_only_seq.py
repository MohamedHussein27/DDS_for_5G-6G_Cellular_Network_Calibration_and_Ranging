import pyuvm
from pyuvm import *
from top_seq_item import *  # Importing your updated randomized item
import cocotb

from top_seq_item import top_item # Ensure this is imported to access the set_fixed_chirp function

class chirp_only_seq(uvm_sequence):
    """
    TC-002: Run 1 frame with a fixed chirp and zero OFDM symbols.
    This acts as the ultimate stress test for the TX Top Datapath.
    """
    def __init__(self, name="chirp_only_seq"):
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
                    
                    req.set_fixed_chirp(f0=30e6, B=200e6)  # Example fixed chirp for frame 1
                    print(f"Step in sequence item", req.FTW_start, req.FTW_step)
                else:
                    req.rst_n      = 1
                    req.dds_enable = 1
                    req.set_fixed_chirp(f0=30e6, B=200e6)  # Example fixed chirp for subsequent frames
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