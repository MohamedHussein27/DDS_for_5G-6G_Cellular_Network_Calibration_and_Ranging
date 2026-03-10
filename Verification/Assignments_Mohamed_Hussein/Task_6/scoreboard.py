from cocotb.queue import *
from cocotb.triggers import *
import cocotb

from transaction import Transaction

class Scoreboard:

    def __init__(self):
        self.mon2sb  = Queue()
        self.passed = 0
        self.failed = 0

    async def run(self):

        while True:

            tr = await self.mon2sb.get()

            expected = self.model(tr)

            if tr.out == expected:
                self.passed += 1
            else:
                self.failed += 1
                print(f"Mismatch: expected={expected} got={tr.out}")

    def model(self, tr):

        if tr.rst_n == 0:
            return 0

        elif tr.op == 0:
            return tr.a + tr.b

        elif tr.op == 1:
            return tr.a ^ tr.b
        
        elif tr.op == 2:
            return tr.a & tr.b

        elif tr.op == 3:
            return tr.a | tr.b


    def report_test_cases(self):

        print("=================================")
        print("TEST REPORT")
        print("PASSED:", self.passed)
        print("FAILED:", self.failed)
        print("=================================")