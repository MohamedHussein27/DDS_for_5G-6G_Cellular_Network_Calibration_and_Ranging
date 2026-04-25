"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: single_band_test.py

    Description:
        This test verifies the TX datapath under Group B (Single-Band Operation) conditions.
        It executes two distinct sequences sequentially to validate features F-01 and F-02:
        1. Constant M_dds (TC-004): Single chirp tone with constant M_dds.
        2. Full M_dds Sweep (TC-005): Full M_dds sweep from 0 to 1023.
        It ensures the system properly generates and maintains single-tone 
        and linearly sweeping chirp bands accurately.
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
from single_band_sequences import single_chirp_tone_seq, full_dds_sweep_seq

@pyuvm.test()
class single_band_test(base_test):
    def build_phase(self):
        super().build_phase()

        # our group B sequences (single-band sequences)
        self.seq_single_tone = single_chirp_tone_seq.create("seq_single_tone")
        self.seq_full_sweep  = full_dds_sweep_seq.create("seq_full_sweep")

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        # 2. Call the universal setup from base_test (Reset)
        await self.run_initial_setup()

        # 3. Explicitly execute the Group B Sequences
        self.logger.info("--- Executing TC-004: Constant M_dds — Single Chirp Tone ---")
        await self.seq_single_tone.start(self.env.agt.sqr)

        self.logger.info("--- Executing TC-005: Full M_dds Sweep — 0 to 1023 ---")
        await self.seq_full_sweep.start(self.env.agt.sqr)

        # 4. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS single tone generation and linear sweeping")
        self.logger.info("---------------------------------------------------------")