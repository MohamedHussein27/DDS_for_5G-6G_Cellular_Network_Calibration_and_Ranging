import cocotb
from cocotb.triggers import RisingEdge, Timer
from transaction import Transaction


class Driver:
    def __init__(self, dut, queue):
        self.dut = dut
        self.queue = queue

    async def run(self):
        for i in range(10):
            txn = Transaction(i, i + 1)

            self.dut.valid.value = 1
            self.dut.data_in.value = txn.data

            await RisingEdge(self.dut.clk)

            self.dut.valid.value = 0

            print(f"[DRIVER] Sent {txn}")

            await self.queue.put(txn)

            await Timer(10, units="ns")