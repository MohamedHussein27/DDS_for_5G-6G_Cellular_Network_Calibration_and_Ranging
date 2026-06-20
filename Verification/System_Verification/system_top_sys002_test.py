"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: sys002_chirp_loopback_test.py

    Description:
        This test verifies the full TX->channel->RX wrapper under SYS-002
        conditions. TX emits an LFM chirp (f0=0MHz, B=200MHz, OFDM disabled)
        looped through an ideal channel. RX processes the chirp via
        FFT -> conjugate multiplier -> IFFT and the resulting range profile
        is checked for a peak at zero delay.
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock import Clock
from pyuvm import *
import logging

import pyuvm

from base_test import base_test
from system_top_sys002_seq import sys002_chirp_loopback_Seq
from system_top_reset_seq import reset_before_frame_seq

@pyuvm.test()
class sys_002_chirp_loopback_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_chirp_loopback = sys002_chirp_loopback_Seq.create("seq_chirp_loopback")
        self.seq_reset_before = reset_before_frame_seq.create("seq_reset_before")

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock
        await self.generate_clock()

        self.logger.info("--- Executing Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.system_agt.sqr)

        # 2. Execute the Chirp Loopback Sequence
        self.logger.info("--- Executing SYS-002: Single-Frame Loopback - Chirp ---")
        await self.seq_chirp_loopback.start(self.env.system_agt.sqr)

        # 3. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : End-to-end chirp loopback - RX range profile peak at zero delay")
        self.logger.info("---------------------------------------------------------")
