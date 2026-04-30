"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
 
    Module: dds_base_test.py
 
    Description:
        Shared base class for all dds block-level tests.
        Inherits from base_test and overrides the agent active/passive
        assignments so only the dds agent is active.
        All individual dds test files import and inherit from dds_base_test.
"""
 
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
import logging
 
from base_test import base_test
 
 
class dds_base_test(base_test):
 
    def build_phase(self):
        super().build_phase()
        # Activate only the dds agent — passivate all others
        ConfigDB().set(self, "env.dds_agt",  "is_active", uvm_active_passive_enum.UVM_ACTIVE)
        ConfigDB().set(self, "env.ifft_agt", "is_active", uvm_active_passive_enum.UVM_PASSIVE)
        ConfigDB().set(self, "env.top_agt",  "is_active", uvm_active_passive_enum.UVM_PASSIVE)
        ConfigDB().set(self, "env.fft_agt",  "is_active", uvm_active_passive_enum.UVM_PASSIVE)
        # Tag so monitors / scoreboards know which DUT block is under test
        ConfigDB().set(self, "*", "VERIF_MODE", "DDS")
 
    # Convenience wrapper used by every run_phase
    async def _run(self, seq):
        await seq.start(self.env.dds_agt.sqr)