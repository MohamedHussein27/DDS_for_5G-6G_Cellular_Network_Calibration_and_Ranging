"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: word_format_sequence.py
"""

import pyuvm
from pyuvm import *
from top_seq_item import top_item
import cocotb
import random

class word_format_normal_seq(uvm_sequence):
    """
    TC-006: Output Word Format — Normal Operation
    Stimulus: Dual-band nominal stimulus (DDS + OFDM active).
    """
    def __init__(self, name="word_format_normal_seq"):
        super().__init__(name)
        
    async def body(self):
        req = top_item.create("req")
        await self.start_item(req)
        
        req.randomize()
        
        # Nominal dual-band operation
        req.rst_n = 1
        req.dds_enable = 1
        req.ofdm_rd_en = 1
        req.cycles = 4096  # Standard FFT window
        
        # Ensure OFDM values are within nominal/safe bounds to avoid intended saturation
        req.ofdm_in_real = random.randint(-16384, 16383) 
        req.ofdm_in_imag = random.randint(-16384, 16383)
        
        await self.finish_item(req)