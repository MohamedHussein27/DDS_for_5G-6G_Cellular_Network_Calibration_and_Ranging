"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ofdm_recovery_test.py

    Description:
        
       
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
from ofdm_recovery_sequences import pure_ofdm_signal_seq, mixed_ofdm_radar_seq, ofdm_sample_count_timing_seq

@pyuvm.test()
class ofdm_recovery_test(base_test):
    def build_phase(self):
        super().build_phase()

        # our group H sequences (regression & stress sequences)
        self.pure_ofdm = pure_ofdm_signal_seq.create("pure_ofdm")
        self.mixed_ofdm_chirp     = mixed_ofdm_radar_seq.create("mixed_ofdm_chirp")
        self.ofdm_sample_count =ofdm_sample_count_timing_seq.create("ofdm_sample_count")

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        # 2. Call the universal setup from base_test (Reset)
        await self.run_initial_setup()

        # 3. Explicitly execute the Group H Sequences
        self.logger.info("--- Executing TC-RX-004: Pure OFDM Recovery After Mixed Frame ---")
        await self.pure_ofdm.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-RX-005: Mixed OFDM+Radar Extraction ---")
        await self.mixed_ofdm_chirp.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-RX-006: OFDM Sample Count & Valid Timing ---")
        await self.ofdm_sample_count.start(self.env.agt.sqr)

        # 4. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Inter-frame contamination, constrained-random stress, guardband limits")
        self.logger.info("---------------------------------------------------------")