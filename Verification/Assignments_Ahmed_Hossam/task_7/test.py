import pyuvm
from pyuvm import *
import cocotb
from env import env
from addition_seq import addition_seq
from xor_seq import xor_seq
from and_seq import and_seq
from or_seq import or_seq

# ----------------------------------------------------------------------------
# BASE TEST
# ----------------------------------------------------------------------------
class base_test(uvm_test):
    def build_phase(self):
        # Build the Environment
        self.env = env.create("env", self)

# ----------------------------------------------------------------------------
# TEST 1: ADDITION
# ----------------------------------------------------------------------------
class test_addition(base_test):
    async def run_phase(self):
        self.raise_objection()
        cocotb.log.info(">>> Running ADDITION Sequence <<<")
        
        seq = addition_seq.create("seq")
        await seq.start(self.env.agent_.sequencer)
        
        self.drop_objection()

# ----------------------------------------------------------------------------
# TEST 2: XOR
# ----------------------------------------------------------------------------
class test_xor(base_test):
    async def run_phase(self):
        self.raise_objection()
        cocotb.log.info(">>> Running XOR Sequence <<<")
        
        seq = xor_seq.create("seq")
        await seq.start(self.env.agent_.sequencer)
        
        self.drop_objection()
 # ----------------------------------------------------------------------------
# TEST 3: XOR
# ----------------------------------------------------------------------------
class test_and(base_test):
    async def run_phase(self):
        self.raise_objection()
        cocotb.log.info(">>> Running AND Sequence <<<")
        
        seq = xor_seq.create("seq")
        await seq.start(self.env.agent_.sequencer)
        
        self.drop_objection()
        
# ----------------------------------------------------------------------------
# TEST 4: OR
# ----------------------------------------------------------------------------
class test_or(base_test):
    async def run_phase(self):
        self.raise_objection()
        cocotb.log.info(">>> Running OR Sequence <<<")
        
        seq = xor_seq.create("seq")
        await seq.start(self.env.agent_.sequencer)
        
        self.drop_objection()



# ----------------------------------------------------------------------------
class test_all(base_test):
    async def run_phase(self):
        self.raise_objection()
        
        cocotb.log.info(">>> Starting ADDITION Sequence <<<")
        seq_add = addition_seq.create("seq_add")
        await seq_add.start(self.env.agent_.sequencer)

        cocotb.log.info(">>> Starting XOR Sequence <<<")
        seq_xor = xor_seq.create("seq_xor")
        await seq_xor.start(self.env.agent_.sequencer)
        cocotb.log.info(">>> Starting and Sequence <<<")
        seq_and = and_seq.create("and_seq")
        await seq_and.start(self.env.agent_.sequencer)
        cocotb.log.info(">>> Starting OR Sequence <<<")
        seq_or = or_seq.create("seq_or")
        await seq_or.start(self.env.agent_.sequencer)
        
        
        self.drop_objection()

# ----------------------------------------------------------------------------
# COCOTB ENTRY POINT
# ----------------------------------------------------------------------------
@cocotb.test()
async def main_test(dut):
   
    ConfigDB().set(None, "*", "DUT", dut)
    
    await uvm_root().run_test("test_all")