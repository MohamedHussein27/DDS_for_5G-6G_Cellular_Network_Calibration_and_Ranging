import cocotb
from cocotb.triggers import FallingEdge


class ALUDriver:
    def __init__(self, dut, seq_queue):
        self.dut = dut
        self.seq_queue = seq_queue

    def initialize(self):
        self.dut.a.value = 0
        self.dut.b.value = 0
        self.dut.op.value = 0
        self.dut.rst_n.value = 1

    async def reset(self, cycles=2):
        cocotb.log.info("Driver: Applying Reset")
        self.dut.rst_n.value = 0
        self.dut.a.value = 0
        self.dut.b.value = 0
        self.dut.op.value = 0

        for _ in range(cycles):
            await FallingEdge(self.dut.clk)

        self.dut.rst_n.value = 1
        cocotb.log.info("Driver: Reset Released")
        await FallingEdge(self.dut.clk)

    async def run(self):
        self.initialize()
        while True:
            tr = await self.seq_queue.get()
            await FallingEdge(self.dut.clk)
            try:
                self.dut.rst_n.value = int(tr.rst_n)
                self.dut.a.value = int(tr.a)
                self.dut.b.value = int(tr.b)
                self.dut.op.value = int(tr.op)
            except (ValueError, TypeError):
                cocotb.log.error(f"Driver: Invalid transaction: {tr}")
                self.initialize()
