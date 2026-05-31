"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: full_sweep_test.py

    Description:
        This test verifies the TX datapath under TC-005 conditions.
        It executes frames with a full M_dds sweep (f0=10MHz, B=200MHz) 
        and randomized OFDM symbols to stress the phase-to-amplitude 
        conversion and datapath mixing logic.
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock import Clock
from pyuvm import *
import logging

import pyuvm

from base_test import base_test
from system_top_normal_seq import normal_Seq
from system_top_reset_seq import reset_before_frame_seq

@pyuvm.test()
class tc_005_full_sweep_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_normal = normal_Seq.create("seq_normal")
        self.seq_reset_before = reset_before_frame_seq.create("seq_reset_before")   

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        self.logger.info("--- Executing Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.system_agt.sqr)

        # 2. Execute the Normal Sweep Sequence
        self.logger.info("--- Executing TC-005: Normal Sweep ---")
        await self.seq_normal.start(self.env.system_agt.sqr)

        # 3. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Normal chirp bandwidth mapping alongside random OFDM data")
        self.logger.info("---------------------------------------------------------")