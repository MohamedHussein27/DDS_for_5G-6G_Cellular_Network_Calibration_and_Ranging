from cocotb.queue import Queue
from agent import ALUAgent
from scoreboard import ALUScoreboard


class ALUEnv:
    def __init__(self):
        # 1. Build Phase: Create the queue
        self.sb_queue = Queue()
        self.agent = None
        self.scoreboard = None

    def connect(self, dut):
        # 2. Connect Phase: Instantiate structural blocks and wire them up
        self.agent = ALUAgent(dut, self.sb_queue)
        self.scoreboard = ALUScoreboard(self.sb_queue)

    def report(self):
        if self.scoreboard:
            self.scoreboard.report()
