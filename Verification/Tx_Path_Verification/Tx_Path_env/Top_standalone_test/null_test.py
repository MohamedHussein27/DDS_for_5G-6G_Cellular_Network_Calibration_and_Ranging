"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: null_test.py

    Description:
        This test verifies the TX datapath under Group A (Null & Boundary) conditions.
        It executes three distinct sequences sequentially to validate features F-01 and F-07:
        1. Null Frame (TC-001): All zeros for both OFDM and DDS paths.
        2. OFDM-Only  (TC-003): Zero chirp, active OFDM transmission.
        3. Chirp-Only (TC-002): Zero OFDM, active DDS chirp sweep.
        It ensures that zeroing one or both inputs produces no artefacts in the active
        band and that valid_out gating functions correctly without data corruption.
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
from null_sequences import null_frame_seq, ofdm_only_seq, chirp_only_seq

@pyuvm.test()
class null_test(base_test):
    def build_phase(self):
        super().build_phase()

        # our group A sequences (null sequences)
        self.seq_null_both  = null_frame_seq.create("seq_null_both")
        self.seq_ofdm_only  = ofdm_only_seq.create("seq_ofdm_only")
        self.seq_chirp_only = chirp_only_seq.create("seq_chirp_only")

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        # 2. Call the universal setup from base_test (Reset)
        await self.run_initial_setup()

        # 3. Explicitly execute the Group A Sequences
        self.logger.info("--- Executing TC-001: Null Frame (Both Zeros) ---")
        await self.seq_null_both.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-003: OFDM-Only (Zero Chirp) ---")
        await self.seq_ofdm_only.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-002: Chirp-Only (Zero OFDM) ---")
        await self.seq_chirp_only.start(self.env.agt.sqr)

        # 4. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : OFDM symbols and chirp null operations")
        self.logger.info("---------------------------------------------------------")