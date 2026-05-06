"""
===============================================================================
DDS Cycle Stress Test
===============================================================================

Description :
    Standalone stress test targeting the DDS cycle counter behavior,
    especially overflow and maximum-count handling.

    This test executes only the cycle stress sequence to simplify
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
# Test 4: Cycle Stress Test
# =============================================================================
@pyuvm.test()
class dds_cyclestress_test(dds_base_test):

    def build_phase(self):
        """Create the cycle stress sequence."""
        
        super().build_phase()

        # Sequence targeting cycle counter stress conditions
        self.cyclestress_seq = dds_cyclestress_seq("cyclestress_seq")
        
    async def run_phase(self):
        """
        Main test execution flow:
            1. Start clock
            2. Apply reset
            3. Run cycle stress sequence
        """

        self.raise_objection()

        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )

        # Start DUT clock
        await self.generate_clock()

        # Apply reset sequence
        await self.run_initial_setup()

        # Start stress sequence
        self.logger.info(
            f"Starting sequence: {self.cyclestress_seq.get_name()}"
        )

        await self.cyclestress_seq.start(self.env.dds_agt.sqr)

        self.logger.info(
            f"finished sequence: {self.cyclestress_seq.get_name()}"
        )

        self.drop_objection()

    def report_phase(self):
        """Print end-of-test summary."""

        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")