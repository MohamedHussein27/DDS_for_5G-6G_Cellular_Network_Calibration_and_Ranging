"""
===============================================================================
DDS Chirp Sweep Standalone Test
===============================================================================

Description :
    This test executes only the DDS chirp sweep sequence in isolation for
    fast debugging and focused functional verification.

    The purpose of this testcase is to:
        1. Validate chirp generation behavior independently.
        2. Reduce simulation runtime during debugging.
        3. Verify DDS sweep functionality without interference from
           additional traffic or mixed scenarios.
        4. Simplify waveform analysis for chirp-related issues.

    The test inherits all common infrastructure from dds_base_test,
    including:
        - Environment construction
        - Clock generation
        - Reset sequencing
        - Agent configuration

===============================================================================
"""

# =============================================================================
# IMPORTS
# =============================================================================

# cocotb core package
import cocotb 
# cocotb synchronization primitives and timing triggers
from cocotb.triggers import * 
# cocotb clock utility
from cocotb.clock import Clock
# cocotb constrained random verification utilities
from cocotb_coverage.crv import *
# pyuvm framework
from pyuvm import * 
# pyuvm decorators/utilities
import pyuvm
# Python logging module
import logging
# =============================================================================
# PROJECT IMPORTS
# =============================================================================

# Base test containing common infrastructure
from dds_base_test import *
# DDS verification sequences
from dds_sequences import * 
# =============================================================================
# TESTCASE: DDS CHIRP SWEEP ONLY TEST
# =============================================================================
@pyuvm.test()
class dds_chirp_only_test(dds_base_test):
    # =========================================================================
    # BUILD PHASE
    # =========================================================================
    def build_phase(self):
        # Build common testbench infrastructure from base class
        super().build_phase()
        # ---------------------------------------------------------------------
        # Create chirp sweep sequence instance
        #
        # This sequence is responsible for generating DDS chirp sweep
        # stimulus transactions.
        # ---------------------------------------------------------------------
        self.chirp_seq = dds_chirpsweep_seq("chirp_seq")
    # =========================================================================
    # RUN PHASE
    # =========================================================================
    async def run_phase(self):
        # Prevent simulation from ending prematurely
        self.raise_objection()
        # ---------------------------------------------------------------------
        # Log test start banner
        # ---------------------------------------------------------------------
        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )
        # ---------------------------------------------------------------------
        # Generate and start DUT clock
        # ---------------------------------------------------------------------
        await self.generate_clock()
        # ---------------------------------------------------------------------
        # Apply DUT reset sequence
        # ---------------------------------------------------------------------
        await self.run_initial_setup()
        # ---------------------------------------------------------------------
        # Log sequence start
        # ---------------------------------------------------------------------
        self.logger.info(
            f"Starting sequence: {self.chirp_seq.get_name()}"
        )
        # ---------------------------------------------------------------------
        # Start chirp sweep sequence on DDS sequencer
        #
        # The sequence sends transactions through the DDS active agent.
        # ---------------------------------------------------------------------
        await self.chirp_seq.start(self.env.dds_agt.sqr)
        # ---------------------------------------------------------------------
        # Log sequence completion
        # ---------------------------------------------------------------------
        self.logger.info(
            f"finished sequence: {self.chirp_seq.get_name()}"
        )
        # ---------------------------------------------------------------------
        # Allow simulation to terminate
        # ---------------------------------------------------------------------
        self.drop_objection() 

    # =========================================================================
    # REPORT PHASE
    # =========================================================================
    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")