"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: normal_chirp_test.py

    Description:
        Test that executes a single frame with the standard linear chirp
        (f0=0, B=200 MHz). No backdoor memory access – the OFDM memory
        is used as pre-loaded. Validates the entire DDS TX datapath.
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
import logging

# Import the base test and the required sequences
from base_test import base_test
from top_seq_item import *
from normal_chirp_seq import normal_chirp_seq
from reset_sequence import reset_before_frame_seq

@pyuvm.test()
class tc_normal_chirp_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_normal_chirp = normal_chirp_seq.create("seq_normal_chirp")
        self.seq_reset_before = reset_before_frame_seq.create("seq_reset_before")

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock (assumed to be defined in base_test)
        await self.generate_clock()

        # 2. Execute the Reset Sequence
        self.logger.info("--- Executing Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.top_agt.sqr)

        # 3. Execute the Normal Chirp Sequence (f0=0, B=200MHz)
        self.logger.info("--- Executing Normal Chirp (Full 200 MHz Bandwidth) Frame ---")
        await self.seq_normal_chirp.start(self.env.top_agt.sqr)

        # 4. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Full-bandwidth chirp (f0=0, B=200MHz),")
        self.logger.info("                 OFDM memory used as-is, no backdoor.")
        self.logger.info("---------------------------------------------------------")