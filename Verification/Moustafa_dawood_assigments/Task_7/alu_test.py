import pyuvm
from pyuvm import *
import cocotb
from alu_env import alu_env

from rst_sequence import rst_sequence
from xor_sequence import xor_sequence
from and_sequence import and_sequence
from or_sequence import or_sequence
from add_sequence import add_sequence
from random_op_pkg import random_op

class base_test(uvm_test):
    def build_phase(self):
        super().build_phase()
        self.env = alu_env("env", self)
        # Pass the DUT to everything
        ConfigDB().set(None, "*", "DUT", cocotb.top)

class rst_test(base_test):
    async def run_phase(self):
        self.raise_objection()
        seq = rst_sequence.create("seq")
        await seq.start(self.env.agt.sqr)
        self.drop_objection()

class add_test(base_test):
    async def run_phase(self):
        self.raise_objection()
        seq = add_sequence.create("seq")
        await seq.start(self.env.agt.sqr)
        self.drop_objection()

class xor_test(base_test):
    async def run_phase(self):
        self.raise_objection()
        seq = xor_sequence.create("seq")
        await seq.start(self.env.agt.sqr)
        self.drop_objection()

class and_test(base_test):
    async def run_phase(self):
        self.raise_objection()
        seq = and_sequence.create("seq")
        await seq.start(self.env.agt.sqr)
        self.drop_objection()

class or_test(base_test):
    async def run_phase(self):
        self.raise_objection()
        seq = or_sequence.create("seq")
        await seq.start(self.env.agt.sqr)
        self.drop_objection()

class random_test(base_test):
    async def run_phase(self):
        self.raise_objection()
        seq = random_op.create("seq")
        await seq.start(self.env.agt.sqr)
        self.drop_objection()

# --- THE MASTER REGRESSION TEST ---
class regression_test(base_test):
    async def run_phase(self):
        self.raise_objection()
        self.logger.info("Starting Grand Regression...")

        seq = rst_sequence.create("rst_seq")
        await seq.start(self.env.agt.sqr)

        seq = add_sequence.create("add_seq")
        await seq.start(self.env.agt.sqr)

        seq = xor_sequence.create("xor_seq")
        await seq.start(self.env.agt.sqr)

        seq = and_sequence.create("and_seq")
        await seq.start(self.env.agt.sqr)

        seq = or_sequence.create("or_seq")
        await seq.start(self.env.agt.sqr)

        seq = random_op.create("rnd_seq")
        await seq.start(self.env.agt.sqr)

        self.logger.info("Grand Regression Complete!")
        self.drop_objection()