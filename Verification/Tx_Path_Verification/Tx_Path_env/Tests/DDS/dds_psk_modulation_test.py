"""
===============================================================================
DDS PSK Modulation Test
===============================================================================

Description :
    Directed verification test targeting PSK modulation functionality
    using the DDS block.

    This test validates:
        - Phase modulation behavior
        - DDS phase transition handling
        - Stable waveform generation during modulation
        - Correct response to PSK sequence updates

    The test executes only the PSK modulation sequence for easier
    debugging and waveform inspection.

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
# Test : DDS PSK Modulation Test
# =============================================================================
@pyuvm.test()
class dds_psk_modulation(dds_base_test):

    def build_phase(self):
        """
        Create the PSK modulation sequence.
        """

        super().build_phase()

        # Sequence generating PSK-modulated DDS output
        self.psk_seq = dds_psk_modulation_seq.create("psk_seq")

    async def run_phase(self):
        """
        Main execution flow:
            1. Start clock
            2. Apply reset
            3. Run PSK modulation sequence
        """

        self.raise_objection()

        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )

        # Start DUT clock
        await self.generate_clock()

        # Apply reset sequence
        await self.run_initial_setup()

        # Start PSK modulation sequence
        self.logger.info(
            f"Starting sequence: {self.psk_seq.get_name()}"
        )

        await self.psk_seq.start(self.env.dds_agt.sqr)

        self.logger.info(
            f"finished sequence: {self.psk_seq.get_name()}"
        )

        self.drop_objection()

    def report_phase(self):
        """Print end-of-test summary."""

        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")