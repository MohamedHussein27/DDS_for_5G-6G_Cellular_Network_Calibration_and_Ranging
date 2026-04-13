from pyuvm import *


class ALU_scoreboard(uvm_component):
    def build_phase(self):
        super().build_phase()
        # FIFO to receive transactions from the agent/monitor
        self.sb_export = uvm_tlm_analysis_fifo("sb_export", self)
        self.passed_count = 0
        self.failed_count = 0

    async def run_phase(self):
        while True:
            # Wait for an item to arrive from the monitor
            item = await self.sb_export.get()

            # Predict the expected output using Python logic
            expected_c, expected_out = self.predict_alu(
                item.a, item.b, item.op)

            # Compare actual vs expected
            if item.out == expected_out and item.c == expected_c:
                self.passed_count += 1
                self.logger.debug(f"PASS: {item}")
            else:
                self.failed_count += 1
                self.logger.error(
                    f"FAIL: {item} | Expected c={expected_c} out={expected_out}")

    def predict_alu(self, a, b, op):
        """Python reference model of the ALU"""
        if op == 0:   # ADD
            res = a + b
        elif op == 1:  # XOR
            res = a ^ b
        elif op == 2:  # AND
            res = a & b
        elif op == 3:  # OR
            res = a | b
        else:
            res = 0

        # Extract carry (bit 4) and 4-bit output (bits 0-3)
        c = (res >> 4) & 1
        out = res & 0xF
        return c, out

    def check_phase(self):
        # Print a summary at the end of the simulation
        if self.failed_count == 0:
            self.logger.info(
                f"TEST PASSED! ({self.passed_count} transactions)")
        else:
            self.logger.error(f"TEST FAILED with {self.failed_count} errors.")
