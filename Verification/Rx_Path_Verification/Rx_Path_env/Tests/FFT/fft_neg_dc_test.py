import cocotb
from pyuvm import *
import pyuvm
from fft_base_test import base_test
from fft_reset_seq import FftResetSeq
from fft_neg_dc_seq import SeqNegDc



import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *

import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *



@pyuvm.test()
class test_fft_neg_dc(base_test):
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        await self.generate_clock()
        await self.run_initial_setup()
        
        self.logger.info("Starting Negative DC Test (-32768)...")
        seq = SeqNegDc("neg_dc_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(5000):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()
