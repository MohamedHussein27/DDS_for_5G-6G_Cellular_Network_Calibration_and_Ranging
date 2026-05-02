import cocotb
from cocotb.triggers import RisingEdge
import pyuvm
from pyuvm import *

from fft_base_test import base_test
from fft_regression_seq import SeqRegression

@pyuvm.test()
class test_fft_master_regression(base_test):
    """
    Test Case: Master Regression Runner
    Description: Executes the regression sequence to test everything in one run.
    """
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        await self.generate_clock()
        await self.run_initial_setup()
        
        # Start the massive regression sequence
        seq = SeqRegression("regression_master")
        await seq.start(self.env.fft_agt.sqr)
        
        # Final safety flush
        self.logger.info("All tests injected. Waiting for final pipeline flush...")
        for _ in range(5000):
            await RisingEdge(cocotb.top.clk)
            
        self.logger.info("Regression Simulation Complete. Dropping objection.")
        self.drop_objection()