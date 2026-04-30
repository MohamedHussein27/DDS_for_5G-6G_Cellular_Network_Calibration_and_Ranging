import cocotb
from pyuvm import *
import pyuvm
from fft_base_test import base_test
from fft_multiframe_seq import fft_two_frame_seq
from fft_reset_seq import FftResetSeq


import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *

import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *




@pyuvm.test()
class test_fft_multiframe(base_test):
    """
    Test Case: Multi-Frame (Two Consecutive Frames)
    Description: Drives two back-to-back 4096-point frames to verify 
    the pipeline can handle continuous data without dropping samples 
    between frames.
    """
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # 1. Start Hardware Clock & Global Reset
        await self.generate_clock()
        await self.run_initial_setup()
        
        # 2. Start the Multi-Frame Sequence
        self.logger.info("Starting Multi-Frame (Continuous) Test...")
        
        # Ensure this matches the class name you used in your sequence file!
        seq = fft_two_frame_seq("two_frame_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # 3. Wait for the pipeline to drain
        # TWO frames = 8192 valid_in cycles. 
        # We must wait >8192 cycles + the RTL pipeline latency.
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(10000):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()