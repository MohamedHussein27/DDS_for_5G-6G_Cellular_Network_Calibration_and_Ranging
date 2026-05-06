"""
===============================================================================
DDS Reset Recovery Test
===============================================================================

Description :
    Directed verification test targeting DDS reset recovery behavior.

    This test validates:
        - Proper recovery after reset assertion
        - Internal state reinitialization
        - Stable operation after reset release
        - DUT robustness during reset transitions

    The test executes only the reset recovery sequence to simplify
    debugging and waveform analysis.

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
# Test 6: DDS Reset Recovery Test
# =============================================================================
@pyuvm.test()
class dds_resetrecovery_only_test(dds_base_test):

    def build_phase(self):
        """
        Create the reset recovery sequence.
        """

        super().build_phase()

        # Sequence targeting reset recovery scenarios
        self.reset_recovery_seq = dds_reset_recovery_seq.create(
            "reset_recovery_seq"
        )

    async def run_phase(self):
        """
        Main execution flow:
            1. Start clock
            2. Apply reset
            3. Run reset recovery sequence
        """

        self.raise_objection()

        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )

        # Start DUT clock
        await self.generate_clock()

        # Apply initial reset sequence
        await self.run_initial_setup()

        # Start reset recovery sequence
        self.logger.info(
            f"Starting sequence: {self.reset_recovery_seq.get_name()}"
        )

        await self.reset_recovery_seq.start(self.env.dds_agt.sqr)

        self.logger.info(
            f"finished sequence: {self.reset_recovery_seq.get_name()}"
        )

        self.drop_objection()

    def report_phase(self):
        """Print end-of-test summary."""

        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")