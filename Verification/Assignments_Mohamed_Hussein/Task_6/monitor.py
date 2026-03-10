from cocotb.queue import *
from cocotb.triggers import *
import cocotb

from transaction import Transaction

class Monitor:

    def __init__(self):
        self.mon2sb  = Queue()
        self.mon2cov = Queue()
        self.tr_mon = Transaction()

    async def run(self, dut):
        await RisingEdge(dut.clk)
        
        while True:

            await RisingEdge(dut.clk)
            await ReadOnly()

            

            self.tr_mon.a      = int(dut.a.value)
            self.tr_mon.b      = int(dut.b.value)
            self.tr_mon.op     = int(dut.op.value)
            self.tr_mon.rst_n  = int(dut.rst_n.value)
            self.tr_mon.out    = int(dut.out.value)
            self.tr_mon.c      = int(dut.c.value)
            cocotb.log.info("[Monitor] begins")
            if self.mon2sb is not None:
                await self.mon2sb.put(self.tr_mon)

            if self.mon2cov is not None:
                await self.mon2cov.put(self.tr_mon)