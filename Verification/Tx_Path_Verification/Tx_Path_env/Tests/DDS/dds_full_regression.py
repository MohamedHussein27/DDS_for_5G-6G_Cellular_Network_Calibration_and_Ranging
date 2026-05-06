"""
===============================================================================
DDS Full Regression Test
===============================================================================

Description :
    Top-level DDS regression test executing multiple directed,
    stress, coverage, and randomized sequences in a single run.

    This test is intended for:
        - Full functional regression
        - Long-run stability validation
        - Cross-feature interaction testing
        - Coverage improvement

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
# Test : Full DDS Regression
# =============================================================================
@pyuvm.test()
class dds_full_test(dds_base_test):

    def build_phase(self):
        """
        Create all regression sequences.
        """

        super().build_phase()

        # Verification mode tag
        ConfigDB().set(self, "*", "VERIF_MODE", "DDS")

        # ---------------------------------------------------------------------
        # Regression sequence collection
        # ---------------------------------------------------------------------

        # Random stress traffic
        self.random_seq = dds_random_seq.create("random_seq")

        # Reset validation sequence
        self.reset_seq = dds_reset.create("reset_seq")

        # Single-tone DDS operation
        self.singletone_seq = dds_singletone_seq.create("singletone_seq")

        # Cycle counter stress sequence
        self.cyclestress_seq = dds_cyclestress_seq.create("cyclestress_seq")

        # Chirp sweep generation
        self.chirpsweep_seq = dds_chirpsweep_seq.create("chirpsweep_seq")

        # FMCW radar scenario
        self.fmcw_seq = dds_fmcw_radar_seq.create("fmcw_seq")

        # PSK modulation test
        self.psk_seq = dds_psk_modulation_seq.create("psk_seq")

        # DC hold operation
        self.dc_hold_seq = dds_dc_hold_seq.create("dc_hold_seq")

        # Maximum functional coverage sequence
        self.max_cov_seq = dds_max_cov_seq.create("max_cov_seq")

        # FFT boundary condition sequence
        self.fft_boundary_seq = dds_fft_boundary_seq.create(
            "dds_fft_boundary_seq"
        )

    async def run_phase(self):
        """
        Main regression execution flow.
        """

        # Ordered regression sequence execution
        test_sequences = [

            # Start with reset
            self.reset_seq,

            # Boundary and directed tests
            self.fft_boundary_seq,
            self.singletone_seq,
            self.cyclestress_seq,

            # Chirp and radar scenarios
            self.chirpsweep_seq,
            self.fmcw_seq,

            # Modulation and coverage tests
            self.psk_seq,
            self.max_cov_seq,
            self.dc_hold_seq,

            # Random stress sequence runs last
            self.random_seq
        ]

        self.raise_objection()

        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )

        # Start DUT clock
        await self.generate_clock()

        # Apply reset sequence
        await self.run_initial_setup()

        # ---------------------------------------------------------------------
        # Execute all regression sequences
        # ---------------------------------------------------------------------
        for seq in test_sequences:

            self.logger.info(
                f"Starting sequence: {seq.get_name()}"
            )

            await seq.start(self.env.dds_agt.sqr)

            self.logger.info(
                f"Finished sequence: {seq.get_name()}"
            )

        # End simulation
        self.drop_objection()

    def report_phase(self):
        """Print regression summary."""

        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")