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
        # 1. Initialize counters for the += operation to work
        self.passed_test_cases = 0
        self.failed_test_cases = 0
        
        self.sc_export = uvm_analysis_export("sc_export", self)
        self.sc_fifo = uvm_tlm_analysis_fifo("sc_fifo", self)
        
    def connect_phase(self):
        self.sc_export.connect(self.sc_fifo.analysis_export)
        
    # 2. Add async here
    async def run_phase(self):
        while True:
            # 3. Add await here
            item = await self.sc_fifo.get()
            
            """************************** TEST CASES **************************"""
            # 4. Check the 'item', not the 'fifo'!
            if item.reset == 0:
                self.reset_test_case(item) # Pass the item to the function
            else:
                match item.op:
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
        self.Golden_out = 0
        self.Golden_c = 0
        if item.out == self.Golden_out and item.c == self.Golden_c:
            self.passed_test_cases += 1
            cocotb.log.info("Reset Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("Reset Test Case Failed ")
    # 5. Receive the 'item' as an argument
    def add_test_case(self, item):
        # Read a and b directly from the item
        full_result = item.a + item.b
        
        self.Golden_c = (full_result >> 4) & 1
        self.Golden_out = full_result & 0xF
   
        # Check against the item's output
        if item.out == self.Golden_out and item.c == self.Golden_c:
            self.passed_test_cases += 1
            cocotb.log.info("Add Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("Add Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{item.out}  GOLDEN C :{self.Golden_c}  YOUR C :{item.c}")
            
    # ... update xor, and, or, and reset functions exactly like add_test_case ...

    def xor_test_case(self,item) :
        self.Golden_out = (item.a ^ item.b) & 0xF
        self.Golden_c   = 0 # Logical operations generate no carry
   
        if item.out == self.Golden_out and item.c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("XOR Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("XOR Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{item.out}  GOLDEN C :{self.Golden_c}  YOUR C :{item.c}")

    def and_test_case(self,item) :
        self.Golden_out = (item.a & item.b) & 0xF
        self.Golden_c   = 0 
   
        if item.out == self.Golden_out and item.c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("AND Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("AND Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{item.out}  GOLDEN C :{self.Golden_c}  YOUR C :{item.c}")

    def or_test_case(self,item) :
        self.Golden_out = (item.a | item.b) & 0xF
        self.Golden_c   = 0 
   
        if item.out == self.Golden_out and item.c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("OR Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("OR Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{item.out}  GOLDEN C :{self.Golden_c}  YOUR C :{item.c}")    