"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: valid_out_gating_test.py

    Description:
        This test verifies the TX datapath under TC-008: valid_out Gating Verification.
        It sends one frame and observes valid_out timing relative to frame_start.
        Verifies valid_out is LOW between frames and during pipeline latency.
        Pass Criteria:
            - valid_out = 1 for exactly 4096 consecutive cycles per frame.
            - valid_out = 0 between frames.
        Pipeline latency L must remain constant (ASS-05 deterministic latency).
        Feature: F-04
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
from valid_out_gating_seq import valid_out_gating_seq
from reset_sequence import reset_before_frame_seq

@pyuvm.test()
class tc_008_valid_out_gating_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_valid_out_gating = valid_out_gating_seq.create("seq_valid_out_gating")
        self.seq_reset_before     = reset_before_frame_seq.create("seq_reset_before")

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock
        await self.generate_clock()

        self.logger.info("--- Executing TC-009: Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.top_agt.sqr)

        # 2. Execute the valid_out Gating Verification Sequence
        self.logger.info("--- Executing TC-008: valid_out Gating Verification ---")
        await self.seq_valid_out_gating.start(self.env.top_agt.sqr)

        # 3. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : valid_out timing relative to frame_start, inter-frame LOW gating")
        self.logger.info(" [PASS CRITERIA]: valid_out count = 4096 AND valid_out = 0 between frames")
        self.logger.info(" [FEATURE]     : F-04 — Pipeline latency L must remain constant (ASS-05 deterministic latency)")
        self.logger.info("---------------------------------------------------------")
