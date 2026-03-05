import cocotb
import cocotb.utils
from cocotb.triggers import RisingEdge, ReadOnly
from transaction import ALUTransaction


class ALUMonitor:
    def __init__(self, dut, sb_queue):
        self.dut = dut
        self.sb_queue = sb_queue

    async def run(self):
        while True:
            await RisingEdge(self.dut.clk)

            # Skip monitoring during reset
            if int(self.dut.rst_n.value) == 0:
                continue

            # Wait until signals settle
            await ReadOnly()

            try:
                if not self.dut.out.value.is_resolvable:
                    cocotb.log.warning(
                        f"Monitor: 'out' is {self.dut.out.value.binstr} at "
                        f"{cocotb.utils.get_sim_time('ns')} ns. Skipping."
                    )
                    continue

                tr = ALUTransaction(
                    int(self.dut.a.value),
                    int(self.dut.b.value),
                    int(self.dut.op.value),
                    int(self.dut.rst_n.value)
                )

                tr.out = int(self.dut.out.value)
                tr.c = int(self.dut.c.value)

                self.sb_queue.put_nowait(tr)

            except ValueError as e:
                cocotb.log.error(f"Monitor caught X/Z value error: {e}")
