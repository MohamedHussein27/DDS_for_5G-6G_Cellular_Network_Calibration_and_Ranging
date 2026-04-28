"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: stress_test.py

    Description:
        This test verifies the TX datapath under Group H (Regression & Stress) conditions.
        It executes three distinct sequences sequentially to validate features F-01, F-02, and F-06:
        1. Back-to-Back 20 Frames (TC-012): Constant transmission with zero gaps.
        2. Long Gap Between Frames (TC-013): Frame 1, massive delay, Frame 2.
        3. 100 Random Frames (TC-014): Full regression with randomized seeds, enforcing 
           OFDM guardband constraints strictly from 200KHz to 210KHz.
"""
import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

# Import the base test and the specific sequences
from base_test import base_test
from stress_sequences import back_to_back_frames_seq, long_gap_frames_seq, random_frames_regression_seq

@pyuvm.test()
class stress_test(base_test):
    def build_phase(self):
        super().build_phase()

        # our group H sequences (regression & stress sequences)
        self.seq_back_to_back = back_to_back_frames_seq.create("seq_back_to_back")
        self.seq_long_gap     = long_gap_frames_seq.create("seq_long_gap")
        self.seq_full_regress = random_frames_regression_seq.create("seq_full_regress")

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        # 2. Call the universal setup from base_test (Reset)
        await self.run_initial_setup()

        # 3. Explicitly execute the Group H Sequences
        self.logger.info("--- Executing TC-012: Back-to-Back 20 Frames — No Gaps ---")
        await self.seq_back_to_back.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-013: Long Gap Between Frames ---")
        await self.seq_long_gap.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-014: 100 Random Frames — Full Regression ---")
        await self.seq_full_regress.start(self.env.agt.sqr)

        # 4. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Inter-frame contamination, constrained-random stress, guardband limits")
        self.logger.info("---------------------------------------------------------")