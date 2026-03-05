import cocotb
from cocotb.queue import Queue
from generator import ALUGenerator
from driver import ALUDriver
from monitor import ALUMonitor


class ALUAgent:
    def __init__(self, dut, sb_queue):
        # 1. Create the internal queue connecting Generator to Driver
        self.gen_queue = Queue()

        # 2. Instantiate the three blocks according to your hierarchy diagram
        self.generator = ALUGenerator(self.gen_queue)
        self.driver = ALUDriver(dut, self.gen_queue)
        self.monitor = ALUMonitor(dut, sb_queue)

    async def run(self):
        # We only start the driver and monitor continuously in the background.
        # The generator will be triggered explicitly by the testbench.
        cocotb.start_soon(self.driver.run())
        cocotb.start_soon(self.monitor.run())
