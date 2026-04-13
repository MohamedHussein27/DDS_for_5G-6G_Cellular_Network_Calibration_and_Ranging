import os
import cocotb
from cocotb.clock import Clock
from cocotb.triggers import Timer
import pyuvm
from pyuvm import *

from ALU_config import ALU_config
from ALU_env import ALU_env

from ALU_add_sequence import ALU_add_sequence
from ALU_and_sequence import ALU_and_sequence
from ALU_or_sequence import ALU_or_sequence
from ALU_xor_sequence import ALU_xor_sequence
from ALU_random_sequence import ALU_random_sequence


# ----------------------------------------------------------------------
# BASE TEST CLASS
# ----------------------------------------------------------------------
class ALU_test(uvm_test):

    def build_phase(self):
        super().build_phase()

        self.env = ALU_env("env", self)
        self.config_test = ALU_config("config_test")

        self.dut = cocotb.top

        self.config_test.vif = self.dut
        ConfigDB().set(self, "*", "CFG", self.config_test)


# ----------------------------------------------------------------------
# SINGLE PYUVM TEST ENTRY POINT
# ----------------------------------------------------------------------
@pyuvm.test()
class ALU_test_runner(ALU_test):

    async def run_phase(self):

        cocotb.start_soon(Clock(self.dut.clk, 10, units="ns").start())

        self.raise_objection()

        test_name = os.getenv("UVM_TESTNAME", "regression")

        # ------------------------------------------------------------
        # ADD TEST
        # ------------------------------------------------------------
        if test_name == "add":
            seq = ALU_add_sequence("seq")
            await seq.start(self.env.ag_env.sr_ag)

        # ------------------------------------------------------------
        # XOR TEST
        # ------------------------------------------------------------
        elif test_name == "xor":
            seq = ALU_xor_sequence("seq")
            await seq.start(self.env.ag_env.sr_ag)

        # ------------------------------------------------------------
        # AND TEST
        # ------------------------------------------------------------
        elif test_name == "and":
            seq = ALU_and_sequence("seq")
            await seq.start(self.env.ag_env.sr_ag)

        # ------------------------------------------------------------
        # OR TEST
        # ------------------------------------------------------------
        elif test_name == "or":
            seq = ALU_or_sequence("seq")
            await seq.start(self.env.ag_env.sr_ag)

        # ------------------------------------------------------------
        # RANDOM TEST
        # ------------------------------------------------------------
        elif test_name == "random":
            seq = ALU_random_sequence("seq")
            await seq.start(self.env.ag_env.sr_ag)

        # ------------------------------------------------------------
        # REGRESSION TEST
        # ------------------------------------------------------------
        elif test_name == "regression":

            regression_plan = [
                ("ADD", ALU_add_sequence),
                ("XOR", ALU_xor_sequence),
                ("AND", ALU_and_sequence),
                ("OR",  ALU_or_sequence),
            ]

            for name, seq_class in regression_plan:

                self.logger.info(f"===== {name} =====")

                seq = seq_class("seq")
                await seq.start(self.env.ag_env.sr_ag)

                await Timer(20, units="ns")

        else:
            self.logger.error(f"Unknown test: {test_name}")
        self.drop_objection()
