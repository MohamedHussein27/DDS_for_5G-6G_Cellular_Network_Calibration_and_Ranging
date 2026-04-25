"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: word_format_test.py

    Description:
        This test verifies the TX datapath under Group C (Word Format & Overflow) conditions.
        It executes three distinct sequences sequentially to validate features F-02 and F-04:
        1. Output Word Format (TC-006): Normal operation, strictly verifying Q5.26 precision.
        2. Overflow / Saturation (TC-007): Worst-case amplitude to check saturation vs wrap-around.
        3. valid_out Gating (TC-008): Verifies pipeline latency and gating control.
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
from word_format_sequences import word_format_normal_seq, overflow_saturation_seq, valid_out_gating_seq

@pyuvm.test()
class word_format_test(base_test):
    def build_phase(self):
        super().build_phase()

        # our group C sequences (word format & overflow sequences)
        self.seq_word_format = word_format_normal_seq.create("seq_word_format")
        self.seq_overflow    = overflow_saturation_seq.create("seq_overflow")
        self.seq_valid_gating= valid_out_gating_seq.create("seq_valid_gating")

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        # 2. Call the universal setup from base_test (Reset)
        await self.run_initial_setup()

        # 3. Explicitly execute the Group C Sequences
        self.logger.info("--- Executing TC-006: Output Word Format — Normal ---")
        await self.seq_word_format.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-007: Overflow / Saturation Behaviour ---")
        await self.seq_overflow.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-008: valid_out Gating Verification ---")
        await self.seq_valid_gating.start(self.env.agt.sqr)

        # 4. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     :  saturation bounds, and valid gating")
        self.logger.info("---------------------------------------------------------")