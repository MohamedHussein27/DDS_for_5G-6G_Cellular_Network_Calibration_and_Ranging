"""
===============================================================================
Sponsor      : Analog Devices, Inc. (ADI)
Institution  : Faculty of Engineering, Ain Shams University
Project      : DDS for 5G/6G Cellular Network Calibration and Ranging

Module       : dds_base_test.py

Description  :
    Shared base test class for all DDS block-level verification tests.

    This class extends the generic base_test and specializes the
    verification environment for standalone DDS verification.

    Main responsibilities:
        - Activate only the DDS agent
        - Configure all remaining agents as passive
        - Set verification mode to DDS
        - Provide a reusable sequence execution helper

===============================================================================
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
import logging

from base_test import base_test


class dds_base_test(base_test):
    """
    DDS-specific base test.

    Used as the parent class for all DDS functional,
    random, stress, and corner-case tests.
    """

    def build_phase(self):
        """
        Configure the environment for DDS-only verification.
        """

        # Build common infrastructure from base_test
        super().build_phase()

        # ---------------------------------------------------------------------
        # Activate only the DDS agent
        # ---------------------------------------------------------------------
        ConfigDB().set(
            self,
            "env.dds_agt",
            "is_active",
            uvm_active_passive_enum.UVM_ACTIVE
        )

        # Remaining agents operate as passive monitors only
        ConfigDB().set(
            self,
            "env.ifft_agt",
            "is_active",
            uvm_active_passive_enum.UVM_PASSIVE
        )

        ConfigDB().set(
            self,
            "env.top_agt",
            "is_active",
            uvm_active_passive_enum.UVM_PASSIVE
        )

        ConfigDB().set(
            self,
            "env.fft_agt",
            "is_active",
            uvm_active_passive_enum.UVM_PASSIVE
        )

        # ---------------------------------------------------------------------
        # Verification mode tag used by monitors/scoreboards
        # ---------------------------------------------------------------------
        ConfigDB().set(self, "*", "VERIF_MODE", "DDS")

    async def _run(self, seq):
        """
        Convenience wrapper for starting DDS sequences.
        """

        await seq.start(self.env.dds_agt.sqr)