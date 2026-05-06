"""
===============================================================================
DDS Reset Test
===============================================================================

Description :
    Directed verification test focused on validating DDS reset behavior.

    This test executes only the reset sequence to verify:
        - Proper DUT initialization
        - Reset recovery behavior
        - Signal clearing after reset assertion
        - Stable operation after reset deassertion

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
# Test 1: DDS Reset Test
# =============================================================================
@pyuvm.test()
class dds_reset_only_test(dds_base_test):

    def build_phase(self):
        """
        Create the reset verification sequence.
        """

        super().build_phase()

        # Sequence targeting reset functionality
        self.reset_seq = dds_reset.create("reset_seq")
    
    async def run_phase(self):
        """
        Main execution flow:
            1. Start clock
            2. Apply reset
            3. Run reset sequence
        """

        self.raise_objection()

        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )

        # Start DUT clock
        await self.generate_clock()

        # Apply initial reset
        await self.run_initial_setup()

        # Start reset sequence
        self.logger.info(
            f"Starting sequence: {self.reset_seq.get_name()}"
        )

        await self.reset_seq.start(self.env.dds_agt.sqr)

        self.logger.info(
            f"Finished sequence: {self.reset_seq.get_name()}"
        )

        self.drop_objection()

    def report_phase(self):
        """Print end-of-test summary."""

        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")