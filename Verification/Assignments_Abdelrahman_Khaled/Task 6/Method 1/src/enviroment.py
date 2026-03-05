from cocotb.queue import Queue
from agent import ALUAgent
from scoreboard import ALUScoreboard


class ALUEnv:
    def __init__(self, dut):
        self.sb_queue = Queue()
        self.agent = ALUAgent(dut, self.sb_queue)
        self.scoreboard = ALUScoreboard(self.sb_queue)

    def start(self):
        self.agent.start()
        self.scoreboard.start()

    def report(self):
        self.scoreboard.report()
