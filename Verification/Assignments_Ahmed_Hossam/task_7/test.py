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
        # Pythonic instantiation
        self.env = env("env", self)

# ----------------------------------------------------------------------------
class test_addition(base_test):
    async def run_phase(self): 
        self.raise_objection() 
        cocotb.log.info(">>> Running ADDITION Sequence <<<")
        seq = addition_seq("seq") 
        await seq.start(self.env.agent_.sqr) 
        self.drop_objection()  

# ----------------------------------------------------------------------------
class test_xor(base_test):
    async def run_phase(self):
        self.raise_objection()
        cocotb.log.info(">>> Running XOR Sequence <<<")
        seq = xor_seq("seq")
        await seq.start(self.env.agent_.sqr)
        self.drop_objection()

# ----------------------------------------------------------------------------
class test_and(base_test):
    async def run_phase(self):
        self.raise_objection()
        cocotb.log.info(">>> Running AND Sequence <<<")
        seq = and_seq("seq")
        await seq.start(self.env.agent_.sqr)
        self.drop_objection()
        
# ----------------------------------------------------------------------------
class test_or(base_test):
    async def run_phase(self):
        self.raise_objection()
        cocotb.log.info(">>> Running OR Sequence <<<")
        seq = or_seq("seq")
        await seq.start(self.env.agent_.sqr)
        self.drop_objection()

# ----------------------------------------------------------------------------
class test_all(base_test):
    async def run_phase(self):
        self.raise_objection()
        
        cocotb.log.info(">>> Starting ADDITION Sequence <<<")
        seq_add = addition_seq("seq_add")
        await seq_add.start(self.env.agent_.sqr)

        cocotb.log.info(">>> Starting XOR Sequence <<<")
        seq_xor = xor_seq("seq_xor")
        await seq_xor.start(self.env.agent_.sqr)
        
        cocotb.log.info(">>> Starting AND Sequence <<<")
        seq_and = and_seq("and_seq")
        await seq_and.start(self.env.agent_.sqr)
        
        cocotb.log.info(">>> Starting OR Sequence <<<")
        seq_or = or_seq("seq_or")
        await seq_or.start(self.env.agent_.sqr)
        
        self.drop_objection()

# ----------------------------------------------------------------------------
@pyuvm.test()
async def main_test(dut):
    # 1. Start the clock in the background (10ns period)
    cocotb.start_soon(Clock(dut.clk, 10, unit="ns").start())
    
   
    
    # 2. Grab the test name and run the UVM phases
    test_name = os.environ.get("UVM_TESTNAME", "test_all")
    await uvm_root().run_test(test_name)
    
    # 3. Once UVM finishes, export the coverage database
    coverage_db.export_to_xml(filename="ALU_coverage.xml")