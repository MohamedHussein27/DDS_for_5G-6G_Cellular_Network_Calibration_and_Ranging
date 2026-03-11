from cocotb.queue import Queue
import cocotb

from driver import ALUDriver
from monitor import ALUMonitor


class ALUAgent:

    def __init__(self, dut):

        self.dut = dut

        self.seq = None
        self.drv = ALUDriver(dut)
        self.mon = ALUMonitor(dut)

        self.s2d = Queue()
        self.m2sb = Queue()
        self.m2cov = Queue()

    def connect(self):

        self.seq.connect(self.s2d)

        self.drv.seq2drv = self.s2d

        self.mon.mon2sb = self.m2sb
        self.mon.mon2cov = self.m2cov

    async def run(self):

        await self.drv.reset()

        cocotb.start_soon(self.seq.run())
        cocotb.start_soon(self.drv.run())
        cocotb.start_soon(self.mon.run())