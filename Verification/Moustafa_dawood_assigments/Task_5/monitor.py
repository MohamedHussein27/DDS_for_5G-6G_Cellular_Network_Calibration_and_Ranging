from cocotb.triggers import RisingEdge, Timer
from transaction import ALUTransaction


class ALUMonitor:

    def __init__(self, dut):

        self.dut = dut
        self.mon2sb = None
        self.mon2cov = None

    async def run(self):

        while True:

            await RisingEdge(self.dut.clk)

            await Timer(1, units="ns")

            item = ALUTransaction()

            item.a = int(self.dut.a.value)
            item.b = int(self.dut.b.value)
            item.op = int(self.dut.op.value)

            c = int(self.dut.c.value)
            out = int(self.dut.out.value)

            item.result = (c << 4) | out

            await self.mon2sb.put(item)
            await self.mon2cov.put(item)