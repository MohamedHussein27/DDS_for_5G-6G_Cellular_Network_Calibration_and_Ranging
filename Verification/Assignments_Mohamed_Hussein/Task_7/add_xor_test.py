import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging


from environment import *
from sequences import *
from test import *


@pyuvm.test()
class add_xor_test(base_test):
    def build_phase(self):
        uvm_factory().set_type_override_by_type(base_sequence, add_xor_sequence)
        super().build_phase()

    def report_phase(self):
        # This will print immediately after the base test finishes running the sequences
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Addition and XOR Operations")
        self.logger.info("---------------------------------------------------------")