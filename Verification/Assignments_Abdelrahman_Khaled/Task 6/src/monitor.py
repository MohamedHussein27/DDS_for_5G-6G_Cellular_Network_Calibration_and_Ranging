import cocotb
from cocotb.triggers import RisingEdge, ReadOnly
from transaction import ALUTransaction


class ALUMonitor:
    def __init__(self, dut, sb_queue):
        self.dut = dut
        self.sb_queue = sb_queue

    async def run(self):
        while True:
            await RisingEdge(self.dut.clk)
            await ReadOnly()

            # Defensive check to ignore uninitialized (X) or high-impedance (Z) states
            if not self.dut.out.value.is_resolvable:
                continue

            tr = ALUTransaction(
                int(self.dut.a.value),
                int(self.dut.b.value),
                int(self.dut.op.value),
                int(self.dut.rst_n.value)
            )
            tr.out = int(self.dut.out.value)
            tr.c = int(self.dut.c.value)

            # Send to scoreboard
            self.sb_queue.put_nowait(tr)
