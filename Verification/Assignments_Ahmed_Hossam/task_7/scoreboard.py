import pyuvm
from pyuvm import *
from seq_item import *
from driver import *
from sequencer import sequencer
from monitor import monitor
from cocotb.triggers import Timer
from cocotb.clock import Clock

class scoreboard(uvm_scoreboard):      
    def build_phase(self):
        self.passed_test_cases = 0
        self.failed_test_cases = 0
        
       
        self.sc_fifo = uvm_tlm_analysis_fifo("sc_fifo", self)
        
    async def run_phase(self):
        while True:
            item = await self.sc_fifo.get()
            
            """************************** TEST CASES **************************"""
            if int(item.reset) == 0:
                self.reset_test_case(item)
            else:
                match int(item.op):
                    case 0:
                        self.add_test_case(item)
                    case 1:
                        self.xor_test_case(item)
                    case 2:
                        self.and_test_case(item)
                    case 3:
                        self.or_test_case(item)
            """******************************************************************"""

    def reset_test_case(self, item):
        val_out = int(item.out)
        val_c = int(item.c)
        
        self.Golden_out = 0
        self.Golden_c = 0
        if val_out == self.Golden_out and val_c == self.Golden_c:
            self.passed_test_cases += 1
            cocotb.log.info("Reset Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("Reset Test Case Failed ")

    def add_test_case(self, item):
        val_a = int(item.a)
        val_b = int(item.b)
        val_out = int(item.out)
        val_c = int(item.c)
        
        full_result = val_a + val_b
        
        self.Golden_c = (full_result >> 4) & 1
        self.Golden_out = full_result & 0xF
   
        # Check against the item's output
        if val_out == self.Golden_out and val_c == self.Golden_c:
            self.passed_test_cases += 1
            cocotb.log.info("Add Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("Add Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{val_out}  GOLDEN C :{self.Golden_c}  YOUR C :{val_c}")
            
    def xor_test_case(self, item) :
        val_a = int(item.a)
        val_b = int(item.b)
        val_out = int(item.out)
        val_c = int(item.c)
        
        self.Golden_out = (val_a ^ val_b) & 0xF
        self.Golden_c   = 0 # Logical operations generate no carry
   
        if val_out == self.Golden_out and val_c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("XOR Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("XOR Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{val_out}  GOLDEN C :{self.Golden_c}  YOUR C :{val_c}")

    def and_test_case(self, item) :
        val_a = int(item.a)
        val_b = int(item.b)
        val_out = int(item.out)
        val_c = int(item.c)
        
        self.Golden_out = (val_a & val_b) & 0xF
        self.Golden_c   = 0 
   
        if val_out == self.Golden_out and val_c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("AND Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("AND Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{val_out}  GOLDEN C :{self.Golden_c}  YOUR C :{val_c}")

    def or_test_case(self, item) :
        val_a = int(item.a)
        val_b = int(item.b)
        val_out = int(item.out)
        val_c = int(item.c)
        
        self.Golden_out = (val_a | val_b) & 0xF
        self.Golden_c   = 0 
   
        if val_out == self.Golden_out and val_c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("OR Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("OR Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{val_out}  GOLDEN C :{self.Golden_c}  YOUR C :{val_c}")    
            
    def report_phase(self):
        cocotb.log.info("========================================")
        cocotb.log.info("         SCOREBOARD FINAL REPORT        ")
        cocotb.log.info("========================================")
        cocotb.log.info(f" PASSED: {self.passed_test_cases}")
        cocotb.log.info(f" FAILED: {self.failed_test_cases}")
        cocotb.log.info("========================================")
        if self.failed_test_cases == 0:
            cocotb.log.info(" TEST STATUS: SUCCESS ")
        else:
            cocotb.log.error(" TEST STATUS: FAILED ")