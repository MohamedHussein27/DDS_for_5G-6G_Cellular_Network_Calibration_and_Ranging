import cocotb
from driver import ALUDriver
from monitor import ALUMonitor
from sequencer import ALUSequencer


class ALUAgent:
    def __init__(self, dut, sb_queue):
        # 1. Instantiate the Sequencer
        self.sequencer = ALUSequencer()

        # 2. Driver pulls from the Sequencer's queue
        self.driver = ALUDriver(dut, self.sequencer.queue)

        # 3. Monitor sends to the Scoreboard queue
        self.monitor = ALUMonitor(dut, sb_queue)

    def start(self):
        cocotb.start_soon(self.driver.run())
        cocotb.start_soon(self.monitor.run())
