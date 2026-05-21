"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: regression_sequence.py
"""

import pyuvm
from pyuvm import *
from top_seq_item import *  # Importing your updated randomized item
import cocotb

class back_to_back_even_frames_seq(uvm_sequence):
    """
    TC-012: Transmit consecutive frames with zero idle cycles between them.
    Modified to 32 frames (even number) per architectural requirements.
    """
    def __init__(self, name="back_to_back_even_frames_seq"):
        super().__init__(name)
        
        self.num_frames = 32
        
    async def body(self):
        for _ in range(self.num_frames):
            req = top_item.create("req")
            await self.start_item(req)
            
            # Generate a new random seed/data set for this frame
            req.randomize() 
            
            # Force valid transmission constraints
            req.rst_n = 1
            req.dds_enable = 1
            req.ofdm_rd_en = 1
            req.cycles = 4096  # Standard 5G FFT Window length
            
            await self.finish_item(req)