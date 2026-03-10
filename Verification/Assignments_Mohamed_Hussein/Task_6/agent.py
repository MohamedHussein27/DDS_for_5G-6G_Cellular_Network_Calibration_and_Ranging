from generator import Generator
from driver import Driver
from monitor import Monitor

from cocotb.queue import Queue
from cocotb.triggers import Event, Timer
import cocotb

class Agent:

    def __init__(self):

        self.gen = Generator()
        self.drv = Driver()
        self.mon = Monitor()

        # will be assigned to the scoreboard and subscriber queues
        self.mon2sb  = Queue()
        self.mon2cov = Queue()

    def connect(self):

        self.drv.gen2drv = self.gen.gen2drv
        self.gen.gen_ack = self.drv.drv_rqt

        self.mon.mon2sb = self.mon2sb
        self.mon.mon2cov = self.mon2cov

    async def run(self, dut, test_type):

        cocotb.start_soon(self.gen.run(test_type))
        cocotb.start_soon(self.drv.run(dut))
        cocotb.start_soon(self.mon.run(dut))
