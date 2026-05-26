"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: single_tone_test.py

    Description:
        This test verifies the TX datapath under TC-004: Single Tone conditions.
        It executes 100 frames with a constant frequency to validate 
        the phase-to-amplitude conversion and Quarter-Wave symmetry logic.
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
from reset_sequence import reset_before_frame_seq
from chirp_only_seq import chirp_only_seq
@pyuvm.test()
class tc_002_chirp_only_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_single_tone = chirp_only_seq.create("seq_single_tone")
        self.seq_reset_before = reset_before_frame_seq.create("seq_reset_before")   

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock (Assuming generate_clock() is defined in base_test)
        await self.generate_clock()

        self.logger.info("--- Executing TC-009: Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.top_agt.sqr)

        # 2. Execute the Single-Tone Sequence
        self.logger.info("--- Executing TC-002: Chirp Only Zero OFDM ---")
        await self.seq_single_tone.start(self.env.top_agt.sqr)

        # 3. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Inter-frame contamination, constrained-random stress, guardband limits")
        self.logger.info("---------------------------------------------------------")