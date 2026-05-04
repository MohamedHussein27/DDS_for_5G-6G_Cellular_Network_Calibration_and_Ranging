"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ofdm_recovery_sequences.py
"""

import pyuvm
from pyuvm import *
from top_seq_item import rx_item
import cocotb
import random
# =========================================================================
# Sequence 1: TC-RX-004 - Pure OFDM Signal Recovery
# =========================================================================
class pure_ofdm_signal_seq(uvm_sequence):
    """
    TC-RX-004: Apply full OFDM signal. Zero chirp reference.
    """
    def _init_(self, name="pure_ofdm_signal_seq"):
        super()._init_(name)
        
    async def body(self):
        # 1. Zero out the Reference RAM
        for i in range(2048):
            req = rx_item.create(f"req_ref_{i}")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.rx_valid_in = 0
            req.ref_wr_en = 1
            req.ref_wr_re = 0
            req.ref_wr_im = 0
            await self.finish_item(req)
            
        # 2. Apply Pure OFDM Data to Channel
        for i in range(4096):
            req = rx_item.create(f"req_rx_{i}")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.ref_wr_en = 0
            req.rx_valid_in = 1
            await self.finish_item(req)


# =========================================================================
# Sequence 2: TC-RX-005 - Mixed OFDM+Radar Extraction
# =========================================================================
class mixed_ofdm_radar_seq(uvm_sequence):
    """
    TC-RX-005: Load valid chirp reference into RAM. Apply combined OFDM+Chirp.
    """
    def _init_(self, name="mixed_ofdm_radar_seq"):
        super()._init_(name)
        
    async def body(self):
        # 1. Load valid randomized chirp reference into RAM
        for i in range(2048):
            req = rx_item.create(f"req_ref_{i}")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.rx_valid_in = 0
            req.ref_wr_en = 1
            await self.finish_item(req)
            
        # 2. Apply combined signal to RX Channel
        for i in range(4096):
            req = rx_item.create(f"req_rx_{i}")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.ref_wr_en = 0
            req.rx_valid_in = 1
            await self.finish_item(req)


# =========================================================================
# Sequence 3: TC-RX-006 - OFDM Sample Count & Valid Timing
# =========================================================================
class ofdm_sample_count_timing_seq(uvm_sequence):
    """
    TC-RX-006: Send one frame to observe ofdm_valid timing and frame gaps.
    """
    def _init_(self, name="ofdm_sample_count_timing_seq"):
        super()._init_(name)
        
    async def body(self):
        # Send exactly one frame (4096 cycles)
        for i in range(4096):
            req = rx_item.create(f"req_rx_{i}")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.ref_wr_en = 0
            req.rx_valid_in = 1
            await self.finish_item(req)
            
        # Drop valid to allow latency/timing observation
        idle = rx_item.create("idle")
        await self.start_item(idle)
        idle.rst_n = 1
        idle.rx_valid_in = 0
        idle.ref_wr_en = 0
        await self.finish_item(idle)