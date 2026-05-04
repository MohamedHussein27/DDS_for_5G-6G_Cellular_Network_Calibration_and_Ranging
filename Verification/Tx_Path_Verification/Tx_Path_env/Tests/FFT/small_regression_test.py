import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *

import pyuvm
# Base Test
from fft_base_test import base_test

# Regression Sequence
from small_regression_seq import SeqSmallRegression


@pyuvm.test()
class test_small_regression(base_test):
    """
    Test Case: Small FFT Regression (3 Sequences)
    Description: Runs a subset of FFT sequences back-to-back.
    """

    def build_phase(self):
        # Configure environment for FFT-only mode
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()

        # 1. Setup
        await self.generate_clock()
        await self.run_initial_setup()

        # 2. Run regression sequence
        self.logger.info("Starting Small FFT Regression...")
        seq = SeqSmallRegression("small_reg_seq")
        await seq.start(self.env.fft_agt.sqr)

        # 3. Pipeline flush (keep for now, but should improve later)
        self.logger.info("Flushing pipeline...")
        for _ in range(5000):
            await RisingEdge(self.dut.clk)

        # 4. End test
        self.logger.info("Test Complete.")
        self.drop_objection()