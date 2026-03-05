import cocotb


class ALUScoreboard:
    def __init__(self, sb_queue):
        self.sb_queue = sb_queue
        self.total = 0
        self.errors = 0

    def start(self):
        cocotb.start_soon(self.run())

    async def run(self):  # FIXED: Was named 'start' previously
        while True:
            tr = await self.sb_queue.get()

            if tr.rst_n == 0:
                continue

            expected = 0
            a, b = tr.a, tr.b

            if tr.op == 0:    # ADD
                expected = (a + b) & 0xF
            elif tr.op == 1:  # XOR
                expected = a ^ b
            elif tr.op == 2:  # AND
                expected = a & b
            elif tr.op == 3:  # OR
                expected = a | b

            if tr.out != expected:
                self.errors += 1
                cocotb.log.error(f"ERROR: Expected={expected} | Got: {tr}")
            else:
                cocotb.log.info(f"PASS: {tr}")

            self.total += 1

    def report(self):
        cocotb.log.info("========== SCOREBOARD REPORT ==========")
        cocotb.log.info(f"Total Checked: {self.total}")
        cocotb.log.info(f"Errors: {self.errors}")
        if self.errors > 0:
            cocotb.log.error("TEST FAILED")
        else:
            cocotb.log.info("TEST PASSED")
