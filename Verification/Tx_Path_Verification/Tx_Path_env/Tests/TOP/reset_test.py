"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: reset_test.py

    Description:
        This test verifies the TX datapath under Group F (Reset & Control) conditions.
        It executes three distinct sequences sequentially to validate feature F-03:
        1. Reset Before Frame Start (TC-009): Validates clean startup from idle.
        2. Reset Mid-Frame (TC-010): Asserts reset during active TX, ensuring immediate halt.
        3. Multiple Consecutive Resets (TC-011): Spam resets to ensure no accumulated state.
        It guarantees the hardware gracefully enters a known zero-state and 
        recovers seamlessly on the next frame without residual corruption.
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
from reset_sequences import reset_before_frame_seq, reset_mid_frame_seq, multiple_resets_seq

@pyuvm.test()
class reset_test(base_test):
    def build_phase(self):
        super().build_phase()

        # our group F sequences (reset & control sequences)
        self.seq_reset_before = reset_before_frame_seq.create("seq_reset_before") 
        self.seq_reset_mid    = reset_mid_frame_seq.create("seq_reset_mid")
        self.seq_reset_multi  = multiple_resets_seq.create("seq_reset_multi")

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        # 2. Call the universal setup from base_test (Reset)
        await self.run_initial_setup() # we can ignore calling this reset task here

        # 3. Explicitly execute the Group F Sequences
        self.logger.info("--- Executing TC-009: Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-010: Reset Mid-Frame ---")
        await self.seq_reset_mid.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-011: Multiple Consecutive Resets ---")
        await self.seq_reset_multi.start(self.env.agt.sqr)

        # 4. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Synchronous hardware recovery and mid-flight halting")
        self.logger.info("---------------------------------------------------------")