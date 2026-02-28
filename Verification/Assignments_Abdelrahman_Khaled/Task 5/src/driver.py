import cocotb
from cocotb.triggers import FallingEdge


class ALUDriver:
    def __init__(self, dut, seq_queue):
        self.dut = dut
        self.seq_queue = seq_queue

    async def run(self):
        while True:
            tr = await self.seq_queue.get()

            # Drive inputs on the FALLING edge so they are stable
            # before the DUT samples them on the next rising edge.
            await FallingEdge(self.dut.clk)

            try:
                self.dut.a.value = int(tr.a)
                self.dut.b.value = int(tr.b)
                self.dut.op.value = int(tr.op)
            except (ValueError, TypeError):
                cocotb.log.error(
                    f"Driver: Received invalid transaction data: {tr}")
                # You might drive a default '0' here to keep simulation running
                self.dut.a.value = 0
                self.dut.b.value = 0
                self.dut.op.value = 0
