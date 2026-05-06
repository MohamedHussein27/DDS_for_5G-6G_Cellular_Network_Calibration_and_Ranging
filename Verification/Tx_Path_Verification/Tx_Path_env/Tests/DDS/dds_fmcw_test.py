"""
===============================================================================
DDS FMCW Radar Test
===============================================================================

Description :
    Functional test targeting an FMCW radar use-case using the DDS block.

    The test executes a realistic radar chirp generation scenario with
    varying sweep parameters to validate DDS behavior in radar-oriented
    applications such as:
        - Frequency sweeping
        - Multi-chirp generation
        - Dynamic FTW updates
        - Continuous waveform generation

===============================================================================
"""

import cocotb 
from cocotb.triggers import * 
from cocotb.clock import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

# Base test and DDS sequences
from dds_base_test import *
from dds_sequences import * 


# =============================================================================
# Test : FMCW Radar Sequence Test
# =============================================================================
@pyuvm.test()
class dds_fmcw_test(dds_base_test):

    def build_phase(self):
        """
        Create the FMCW radar sequence.
        """

        super().build_phase()

        # Sequence generating FMCW radar chirps
        self.fmcw_seq = dds_fmcw_radar_seq.create("fmcw_seq")

    async def run_phase(self):
        """
        Main execution flow:
            1. Start clock
            2. Apply reset
            3. Run FMCW radar sequence
        """

        self.raise_objection()

        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )

        # Start DUT clock
        await self.generate_clock()

        # Apply reset sequence
        await self.run_initial_setup()

        # Start FMCW sequence
        self.logger.info(
            f"Starting sequence: {self.fmcw_seq.get_name()}"
        )

        await self.fmcw_seq.start(self.env.dds_agt.sqr)

        self.logger.info(
            f"finished sequence: {self.fmcw_seq.get_name()}"
        )

        self.drop_objection()

    def report_phase(self):
        """Print end-of-test summary."""

        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")