import pyuvm
import os
from pyuvm import *
import cocotb
from env import env
from addition_seq import addition_seq
from xor_seq import xor_seq
from and_seq import and_seq
from or_seq import or_seq
from cocotb_coverage.coverage import coverage_db
from cocotb.triggers import Timer
from cocotb.clock import Clock

# ----------------------------------------------------------------------------
class base_test(uvm_test):
    def build_phase(self):
        self.env = env.create("env", self)

# ----------------------------------------------------------------------------
class test_addition(base_test):
    async def run_phase(self): # NO 'phase' argument!
        self.raise_objection() # Use self.raise_objection
        cocotb.log.info(">>> Running ADDITION Sequence <<<")
        seq = addition_seq.create("seq")
        await seq.start(self.env.agent_.sequencer)
        self.drop_objection()  # Use self.drop_objection

# ----------------------------------------------------------------------------
class test_xor(base_test):
    async def run_phase(self):
        self.raise_objection()
        cocotb.log.info(">>> Running XOR Sequence <<<")
        seq = xor_seq.create("seq")
        await seq.start(self.env.agent_.sequencer)
        self.drop_objection()

# ----------------------------------------------------------------------------
class test_and(base_test):
    async def run_phase(self):
        self.raise_objection()
        cocotb.log.info(">>> Running AND Sequence <<<")
        seq = and_seq.create("seq")
        await seq.start(self.env.agent_.sequencer)
        self.drop_objection()
        
# ----------------------------------------------------------------------------
class test_or(base_test):
    async def run_phase(self):
        self.raise_objection()
        cocotb.log.info(">>> Running OR Sequence <<<")
        seq = or_seq.create("seq")
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
        
        cocotb.log.info(">>> Starting AND Sequence <<<")
        seq_and = and_seq.create("and_seq")
        await seq_and.start(self.env.agent_.sequencer)
        
        cocotb.log.info(">>> Starting OR Sequence <<<")
        seq_or = or_seq.create("seq_or")
        await seq_or.start(self.env.agent_.sequencer)
        
        self.drop_objection()

# ----------------------------------------------------------------------------
@cocotb.test()
async def main_test(dut):
    # 1. Start the clock in the background (10ns period)
    # Your DUT relies on this to update out_data!
    cocotb.start_soon(Clock(dut.clk, 10, units="ns").start())
    
    # 2. Pass the DUT to the UVM database
    ConfigDB().set(None, "*", "DUT", dut)
    
    # 3. Grab the test name and run the UVM phases
    test_name = os.environ.get("UVM_TESTNAME", "test_all")
    await uvm_root().run_test(test_name)
    
    # 4. Once UVM finishes, export the coverage database
    coverage_db.export_to_xml(filename="ALU_coverage.xml")