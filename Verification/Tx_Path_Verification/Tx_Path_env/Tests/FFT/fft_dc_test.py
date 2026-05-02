import cocotb
from pyuvm import *
import pyuvm
from fft_base_test import base_test
from fft_reset_seq import FftResetSeq
from fft_dc_seq import SeqDc



import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *

import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *

@pyuvm.test()
class test_fft_dc(base_test):
    """
    Test Case: DC Input Test
    Description: Drives a constant signal (29490) across all 4096 samples 
    to verify the DC bin (Bin 0) accumulates correctly and no other 
    frequency bins show phantom energy.
    """
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # 1. Start Hardware Clock & Global Reset
        await self.generate_clock()
        await self.run_initial_setup()
        
        # 2. Start the DC Sequence
        self.logger.info("Starting DC Input Test (Constant 29490)...")
        
        # Ensure this matches the class name in your sequence file
        seq = SeqDc("dc_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # 3. Wait for the pipeline to flush
        # 4096 cycles for data to enter + RTL pipeline latency
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(5000):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()
