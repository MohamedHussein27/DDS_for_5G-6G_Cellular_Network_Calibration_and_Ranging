"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: radar_sequences.py
"""

import pyuvm
from pyuvm import *
from top_seq_item import rx_item
import cocotb
import random
# =========================================================================
# Sequence 4: TC-RX-007 - Radar Range Profile Loaded Reference
# =========================================================================
class radar_range_profile_seq(uvm_sequence):
    """
    TC-RX-007: Pre-load 2048-entry RAM. Apply matching chirp to channel.
    """
    def _init_(self, name="radar_range_profile_seq"):
        super()._init_(name)
        
    async def body(self):
        ref_re = [random.randint(-10000, 10000) for _ in range(2048)]
        ref_im = [random.randint(-10000, 10000) for _ in range(2048)]
        
        # 1. Load Reference RAM
        for i in range(2048):
            req = rx_item.create(f"req_ref_{i}")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.rx_valid_in = 0
            req.ref_wr_en = 1
            req.ref_wr_re = ref_re[i]
            req.ref_wr_im = ref_im[i]
            await self.finish_item(req)
            
        # 2. Stream Channel (Inject matching chirp into the upper 2048 bins)
        for i in range(4096):
            req = rx_item.create(f"req_rx_{i}")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.ref_wr_en = 0
            req.rx_valid_in = 1
            
            # Simulate perfectly matched radar target in the upper band
            if i >= 2048:
                req.rx_in_re = ref_re[i - 2048]
                req.rx_in_im = ref_im[i - 2048]
                
            await self.finish_item(req)


# =========================================================================
# Sequence 5: TC-RX-008 - Radar Sample Count per Frame
# =========================================================================
class radar_sample_count_seq(uvm_sequence):
    """
    TC-RX-008: Send one frame with valid chirp to verify radar_valid timing.
    """
    def _init_(self, name="radar_sample_count_seq"):
        super()._init_(name)
        
    async def body(self):
        # Load a dummy chirp to activate multiplier
        for i in range(2048):
            req = rx_item.create(f"req_ref_{i}")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.rx_valid_in = 0
            req.ref_wr_en = 1
            await self.finish_item(req)

        # Send exactly one frame
        for i in range(4096):
            req = rx_item.create(f"req_rx_{i}")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.ref_wr_en = 0
            req.rx_valid_in = 1
            await self.finish_item(req)
            
        # Drop valid for timing measurement
        idle = rx_item.create("idle")
        await self.start_item(idle)
        idle.rst_n = 1
        idle.rx_valid_in = 0
        idle.ref_wr_en = 0
        await self.finish_item(idle)