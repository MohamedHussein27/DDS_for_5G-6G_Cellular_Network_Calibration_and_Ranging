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

            # Wait for signals to settle in the current clock cycle
            await ReadOnly()

            try:
                # Use .is_resolvable to check status before casting
                if not self.dut.out.value.is_resolvable:
                    cocotb.log.warning(
                        f"Monitor: 'out' is {self.dut.out.value.binstr} at {cocotb.utils.get_sim_time('ns')} ns. Skipping.")
                    continue

                tr = ALUTransaction(
                    int(self.dut.a.value),
                    int(self.dut.b.value),
                    int(self.dut.op.value)
                )
                tr.out = int(self.dut.out.value)
                tr.c = int(self.dut.c.value)

                self.sb_queue.put_nowait(tr)

            except ValueError as e:
                # Fallback: catch cases where cast fails
                cocotb.log.error(f"Monitor caught X/Z value error: {e}")
