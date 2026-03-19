import pyuvm
from pyuvm import *
from alu_pkg import opcode_t

class alu_scoreboard(uvm_scoreboard):
    def build_phase(self):
        super().build_phase()
        self.sb_export = uvm_analysis_export("sb_export", self)
        self.sb_fifo = uvm_tlm_analysis_fifo("sb_fifo", self)
        
        self.error_count = 0
        self.correct_count = 0
        self.unique_errors = set() 

    def connect_phase(self):
        super().connect_phase()
        self.sb_export = self.sb_fifo.analysis_export

    async def run_phase(self):
        while True:
            item = await self.sb_fifo.get() 
            
            if not item.rst_n:
                expected = 0
            else:
                if item.op == opcode_t.ADD_op:   expected = (item.a + item.b)
                elif item.op == opcode_t.XOR_op: expected = item.a ^ item.b
                elif item.op == opcode_t.AND_op: expected = item.a & item.b
                elif item.op == opcode_t.OR_op:  expected = item.a | item.b
                else: expected = 0

            if item.result != expected:
                self.error_count += 1
                error_sig = f"A:{item.a:2} B:{item.b:2} Op:{opcode_t(item.op).name:6} | Exp:{expected:2} Got:{item.result:2}"
                self.unique_errors.add(error_sig)
            else:
                self.correct_count += 1

    def report_phase(self):
        self.logger.info("==================================================")
        self.logger.info("       ALU GRAND REGRESSION SCOREBOARD REPORT     ")
        self.logger.info("==================================================")
        self.logger.info(f"Total Correct Transactions: {self.correct_count}")
        self.logger.info(f"Total Gross Errors:         {self.error_count}")
        self.logger.info(f"Total UNIQUE Errors:        {len(self.unique_errors)}")
        self.logger.info("--------------------------------------------------")
        
        if len(self.unique_errors) > 0:
            self.logger.error("--- UNIQUE ERROR SUMMARY ---")
            for err in self.unique_errors:
                self.logger.error(err)
            self.logger.error("----------------------------")
            ConfigDB().set(None, "", "TEST_RESULT", "FAIL")
        else:
            self.logger.info("TEST PASSED WITH ZERO ERRORS!")

    def check_phase(self):
        
        if len(self.unique_errors) > 0:
            self.logger.error(f"TEST COMPLETED: Scoreboard successfully caught {len(self.unique_errors)} injected faults.")