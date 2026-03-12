from transaction import * 
from cocotb.triggers import * 
import cocotb
import cocotb.queue

class scoreboard () :
    def __init__(self ,name = "SCOREBOARD"): 
       self.name                = name
       self.t_score             = transaction()
       self.score_mail          = cocotb.queue.Queue()
       self.passed_test_cases   = 0
       self.failed_test_cases   = 0
       self.Golden_c     = 0
       self.Golden_out   = 0

    async def run_scoreboard (self) : 
        while(True):
            self.t_score = transaction()
            self.t_score = await self.score_mail.get() 
            cocotb.log.info("[Scoreboard] receiving from monitor..... ") 
            
            """************************** TEST CASES **************************"""
            if self.t_score.reset == 0:
                self.reset_test_case()
            elif self.t_score.op == 0 and self.t_score.reset == 1:
                self.add_test_case()    
            elif self.t_score.op == 1 and self.t_score.reset == 1:
                self.xor_test_case()
            elif self.t_score.op == 2 and self.t_score.reset == 1:
                self.and_test_case()
            elif self.t_score.op == 3 and self.t_score.reset == 1:
                self.or_test_case()
            """******************************************************************"""

    def reset_test_case (self) :
        self.Golden_out = 0
        self.Golden_c   = 0
        if self.t_score.out == self.Golden_out and self.t_score.c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("Reset Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.info("Reset Test Case Failed ")

    def add_test_case(self) :
        # Calculate full 5-bit result
        full_result = self.t_score.a + self.t_score.b
        
        # Shift right by 4 to grab the 5th bit (the carry)
        self.Golden_c = (full_result >> 4) & 1
        
        # Bitwise AND with 0xF (1111 in binary) to isolate the bottom 4 bits
        self.Golden_out = full_result & 0xF
   
        if self.t_score.out == self.Golden_out and self.t_score.c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("Add Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("Add Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{self.t_score.out}  GOLDEN C :{self.Golden_c}  YOUR C :{self.t_score.c}")

    def xor_test_case(self) :
        self.Golden_out = (self.t_score.a ^ self.t_score.b) & 0xF
        self.Golden_c   = 0 # Logical operations generate no carry
   
        if self.t_score.out == self.Golden_out and self.t_score.c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("XOR Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("XOR Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{self.t_score.out}  GOLDEN C :{self.Golden_c}  YOUR C :{self.t_score.c}")

    def and_test_case(self) :
        self.Golden_out = (self.t_score.a & self.t_score.b) & 0xF
        self.Golden_c   = 0 
   
        if self.t_score.out == self.Golden_out and self.t_score.c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("AND Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("AND Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{self.t_score.out}  GOLDEN C :{self.Golden_c}  YOUR C :{self.t_score.c}")

    def or_test_case(self) :
        self.Golden_out = (self.t_score.a | self.t_score.b) & 0xF
        self.Golden_c   = 0 
   
        if self.t_score.out == self.Golden_out and self.t_score.c == self.Golden_c :
            self.passed_test_cases += 1
            cocotb.log.info("OR Test Case Passed ")
        else:
            self.failed_test_cases += 1
            cocotb.log.error("OR Test Case Failed ")
            cocotb.log.info(f"GOLDEN OUT :{self.Golden_out}  YOUR OUT :{self.t_score.out}  GOLDEN C :{self.Golden_c}  YOUR C :{self.t_score.c}")