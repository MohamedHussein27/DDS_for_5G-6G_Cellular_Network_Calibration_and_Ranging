"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: OFDM_Only_test.py

    Description:
        This test verifies the TX datapath under TC-003: OFDM-Only conditions.
        It executes frames with zero chirp (M_dds = 0) and randomized OFDM
        symbols to validate the datapath mapping and IFFT operations.
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
from OFDM_only_seq import ofdm_only_seq
from reset_sequence import reset_before_frame_seq

@pyuvm.test()
class tc_003_ofdm_only_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_ofdm_only = ofdm_only_seq.create("seq_ofdm_only")
        self.seq_reset_before = reset_before_frame_seq.create("seq_reset_before")   

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock (Assuming generate_clock() is defined in base_test)
        await self.generate_clock()

        self.logger.info("--- Executing TC-009: Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.top_agt.sqr)

        # 2. Execute the OFDM-Only Sequence
        self.logger.info("--- Executing TC-003: OFDM-Only (Zero Chirp) Frames ---")
        await self.seq_ofdm_only.start(self.env.top_agt.sqr)

        # 3. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : OFDM random symbol mapping, datapath routing, zero-chirp math")
        self.logger.info("---------------------------------------------------------")