"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: long_gap_between_frames_test.py

    Description:
        This test verifies the TX datapath under TC-013: Long Gap Between Frames.
        It sends frame 1, inserts 10,000 idle cycles (no frame_start), then sends
        frame 2 with a completely different stimulus.
        Expected: Frame 2 contains no stale pipeline state from frame 1.
        PASS: Frame 2 outputs exactly match the Python reference (Error = 0).
        Priority: MED
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
from top_seq_item import *
from long_gap_between_frames_seq import long_gap_between_frames_seq
from reset_sequence import reset_before_frame_seq

@pyuvm.test()
class tc_013_long_gap_between_frames_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_long_gap     = long_gap_between_frames_seq.create("seq_long_gap")
        self.seq_reset_before = reset_before_frame_seq.create("seq_reset_before")

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock
        await self.generate_clock()

        self.logger.info("--- Executing TC-009: Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.top_agt.sqr)

        # 2. Execute the Long Gap Between Frames Sequence
        self.logger.info("--- Executing TC-013: Long Gap Between Frames (10,000 idle cycles) ---")
        await self.seq_long_gap.start(self.env.top_agt.sqr)

        # 3. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Stale pipeline state after 10,000-cycle idle gap between frames")
        self.logger.info(" [PASS CRITERIA]: Frame 2 outputs exactly match the Python reference (Error = 0)")
        self.logger.info(" [PRIORITY]    : MED")
        self.logger.info("---------------------------------------------------------")
