import cocotb
from pyuvm import *
import pyuvm
from fft_base_test import base_test
from fft_reset_seq import FftResetSeq
from fft_ten_sample_seq import fft_ten_sample_seq


import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *

import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *



@pyuvm.test()
class test_fft_ten_samples(base_test):
    """
    Test Case: Ten Sample Input with Zero-Padding
    Description: Drives 10 custom real-valued samples followed by zeros 
    to verify the FFT's transient response and frequency spreading.
    """
    def build_phase(self):
        # Set the mode so the Environment knows we are verifying the FFT block
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # 1. Start Hardware Clock & Global Reset from base_test
        await self.generate_clock()
        await self.run_initial_setup()
        
        # 2. Start the Ten-Sample Sequence
        self.logger.info("Starting Ten-Sample Input (Zero-Padded) Test...")
        
        # Ensure the sequence name matches your sequence file class name
        seq = fft_ten_sample_seq("ten_sample_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # 3. Wait for the pipeline to flush
        # For a 4096-point FFT, we wait >4096 cycles to ensure the 
        # result for our 10 samples propagates through all internal stages.
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(4500):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()