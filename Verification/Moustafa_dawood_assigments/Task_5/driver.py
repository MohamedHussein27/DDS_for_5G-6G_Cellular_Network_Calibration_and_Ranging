from cocotb.triggers import FallingEdge, RisingEdge
import cocotb
import logging

class ALUDriver:

    def __init__(self, dut, queue):
        self.dut = dut
        self.queue = queue

    async def reset(self):
        logging.info("DRIVER Asserting Reset")
        self.dut.rst_n.value = 0

        for _ in range(2):
            await RisingEdge(self.dut.clk)

        self.dut.rst_n.value = 1
        logging.info("DRIVER Reset Done")
        await RisingEdge(self.dut.clk)

    async def run(self):

        while True:

            item = await self.queue.get()

            await FallingEdge(self.dut.clk)

            self.dut.a.value = item.a
            self.dut.b.value = item.b
            self.dut.op.value = item.op
