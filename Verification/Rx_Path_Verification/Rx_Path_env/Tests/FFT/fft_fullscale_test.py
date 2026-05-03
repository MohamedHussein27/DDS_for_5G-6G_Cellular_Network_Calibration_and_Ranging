import cocotb
from pyuvm import *
import pyuvm
from fft_base_test import base_test
from fft_reset_seq import FftResetSeq
from fft_fullscale_seq import SeqFullScale



import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *

import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *


@pyuvm.test()
class test_fft_full_scale(base_test):
    """
    Test Case: Full-Scale Nyquist Stress Test
    Description: Pushes alternating max/min values (+32767, -32768) 
    to test the highest possible frequency (Nyquist) and verify 
    the internal multipliers and adders do not overflow.
    """
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        await self.generate_clock()
        await self.run_initial_setup()
        
        self.logger.info("Starting Full-Scale Nyquist Test...")
        
        # Make sure this matches the class name in your sequence file!
        seq = SeqFullScale("full_scale_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(5000):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()
