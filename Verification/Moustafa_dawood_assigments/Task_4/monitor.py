from cocotb.triggers import RisingEdge


class Monitor:
    def __init__(self, dut, out_queue):
        self.dut = dut
        self.out_queue = out_queue

    async def run(self):
        while True:
            await RisingEdge(self.dut.clk)

            if self.dut.ready.value == 1:
                data = int(self.dut.data_out.value)
                print(f"[MONITOR] Captured DATA={data}")
                await self.out_queue.put(data)