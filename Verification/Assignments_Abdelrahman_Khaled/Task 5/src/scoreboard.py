import logging
import cocotb


class ALUScoreboard:
    def __init__(self, sb_queue):
        self.sb_queue = sb_queue
        self.total = 0
        self.errors = 0

    def start(self):
        cocotb.start_soon(self.run())

    async def run(self):
        while True:
            tr = await self.sb_queue.get()

            expected = 0
            carry = 0
            a, b = tr.a, tr.b

            if tr.op == 0:  # ADD
                result = a + b
                expected = result & 0xF
                carry = (result >> 4) & 1
            elif tr.op == 1:  # XOR
                expected = a ^ b
            elif tr.op == 2:  # AND
                expected = a & b
            elif tr.op == 3:  # OR
                expected = a | b

            # Better logging using cocotb.log
            if tr.out != expected:
                self.errors += 1
                cocotb.log.error(f"ERROR: {tr} Expected={expected}")
            else:
                cocotb.log.info(f"PASS: {tr}")

            self.total += 1

    def report(self):
        cocotb.log.info("========== SCOREBOARD REPORT ==========")
        cocotb.log.info(f"Total Transactions: {self.total}")
        cocotb.log.info(f"Errors: {self.errors}")

        if self.errors > 0:
            cocotb.log.error("TEST FAILED")
        else:
            cocotb.log.info("TEST PASSED")
