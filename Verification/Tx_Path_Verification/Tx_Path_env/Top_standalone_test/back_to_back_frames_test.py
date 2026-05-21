"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: back_to_back_frames_test.py

    Description:
        This test verifies the TX datapath under TC-012: Back-to-Back 20 Frames — No Gaps.
        It transmits 20 consecutive frames with zero idle cycles between them.
        Each frame uses a different random OFDM seed to stress inter-frame isolation.
        Expected: No inter-frame contamination observable at tx_out for any of the 20 frames.
        PASS: 0 errors across all 20 frames compared to the fixed-point Python model.
        Priority: HIGH
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
from back_to_back_frames_seq import back_to_back_frames_seq
from reset_sequence import reset_before_frame_seq

@pyuvm.test()
class tc_012_back_to_back_frames_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_back_to_back = back_to_back_frames_seq.create("seq_back_to_back")
        self.seq_reset_before = reset_before_frame_seq.create("seq_reset_before")

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock
        await self.generate_clock()

        self.logger.info("--- Executing TC-009: Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.top_agt.sqr)

        # 2. Execute the Back-to-Back 20 Frames Sequence
        self.logger.info("--- Executing TC-012: Back-to-Back 20 Frames — No Gaps ---")
        await self.seq_back_to_back.start(self.env.top_agt.sqr)

        # 3. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Inter-frame contamination, back-to-back stress with unique OFDM seeds")
        self.logger.info(" [PASS CRITERIA]: 0 errors across all 20 frames compared to fixed-point Python model")
        self.logger.info(" [PRIORITY]    : HIGH")
        self.logger.info("---------------------------------------------------------")
