"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: 

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
from radar_sequences import  radar_range_profile_seq, radar_sample_count_seq

@pyuvm.test()
class ofdm_recovery_test(base_test):
    def build_phase(self):
        super().build_phase()

        # our group H sequences (regression & stress sequences)
        self.radar_range = radar_range_profile_seq.create("radar_range")
        self.radar_sample = radar_sample_count_seq.create("radar_sample")
        

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        # 2. Call the universal setup from base_test (Reset)
        await self.run_initial_setup()

        # 3. Explicitly execute the Group H Sequences
        self.logger.info("--- Executing TC-RX-007: Radar Range Profile - Loaded Reference ---")
        await self.radar_range.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-RX-008: Radar Sample Count per Frame ---")
        await self.radar_sample.start(self.env.agt.sqr)

        

        # 4. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Inter-frame contamination, constrained-random stress, guardband limits")
        self.logger.info("---------------------------------------------------------")