from cocotb.triggers import FallingEdge


class ALUDriver:

    def __init__(self, dut):
        self.dut = dut
        self.seq2drv = None

    async def reset(self):

        self.dut.rst_n.value = 0

        for _ in range(2):
            await FallingEdge(self.dut.clk)

        self.dut.rst_n.value = 1

    async def run(self):

        while True:

            item = await self.seq2drv.get()

            await FallingEdge(self.dut.clk)

            self.dut.a.value = item.a
            self.dut.b.value = item.b
            self.dut.op.value = item.op