"""
===============================================================================
DDS Random Test
===============================================================================

Description :
    Constrained-random DDS verification test.

    This test executes randomized DDS transactions to explore a wide
    range of operating scenarios and uncover unexpected corner cases.

    The test helps validate:
        - General DUT robustness
        - Random FTW combinations
        - Random cycle configurations
        - Stable operation under varying inputs

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
# Test 2: DDS Random Test
# =============================================================================
@pyuvm.test()
class dds_random(dds_base_test):

    def build_phase(self):
        """
        Create the random stimulus sequence.
        """

        super().build_phase()

        # Randomized DDS stimulus sequence
        self.random_seq = dds_random_seq("random_seq")
        
    async def run_phase(self):
        """
        Main execution flow:
            1. Start clock
            2. Apply reset
            3. Run random sequence
        """

        self.raise_objection()

        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )

        # Start DUT clock
        await self.generate_clock()

        # Apply reset sequence
        await self.run_initial_setup()

        # Start random sequence
        self.logger.info(
            f"Starting sequence: {self.random_seq.get_name()}"
        )

        await self.random_seq.start(self.env.dds_agt.sqr)

        self.logger.info(
            f"Finished sequence: {self.random_seq.get_name()}"
        )

        self.drop_objection()