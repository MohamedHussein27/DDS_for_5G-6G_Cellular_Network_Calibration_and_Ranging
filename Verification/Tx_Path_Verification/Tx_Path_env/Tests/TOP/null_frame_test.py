"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: null_frame_test.py

    Description:
        This test verifies the TX datapath under TC-001: Null Frame conditions.
        It executes frames with zero chirp (M_dds = 0) and an entirely zeroed-out
        OFDM RAM to prove there are no DC offsets, rounding errors, or uninitialized 
        memory leaks in the datapath. The expected output is a perfect array of zeros
        for exactly 4096 valid clock cycles.
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
import logging

# Import the base test and the specific sequences
from base_test import base_test
from null_frame_seq import null_frame_seq
from reset_sequence import reset_before_frame_seq

@pyuvm.test()
class tc_001_null_frame_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_null_frame = null_frame_seq.create("seq_null_frame")
        self.seq_reset_before = reset_before_frame_seq.create("seq_reset_before")   

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        self.logger.info("--- Executing TC-009: Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.top_agt.sqr)

        # 2. Execute the Null Frame Sequence
        self.logger.info("--- Executing TC-001: Null Frame (All Zeros) ---")
        await self.seq_null_frame.start(self.env.top_agt.sqr)

        # 3. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Absolute zero datapath integrity, 4096 valid_out count")
        self.logger.info("---------------------------------------------------------")