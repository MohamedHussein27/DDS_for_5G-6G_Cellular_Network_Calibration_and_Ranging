import cocotb
from pyuvm import *
import pyuvm
from fft_base_test import base_test
from fft_impulse_seq import SeqImpulse
from fft_reset_seq import FftResetSeq

import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *

import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *



@pyuvm.test()
class test_fft_impulse(base_test):
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # Start Hardware Clock & Reset
        await self.generate_clock()
        await self.run_initial_setup()
        
        # Start the Impulse Sequence
        self.logger.info("Starting Impulse Response (Delta) Test...")
        seq = SeqImpulse("impulse_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # Wait for the pipeline to flush
        # This is critical to see the frequency response results in the scoreboard
        for _ in range(4500):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.drop_objection()