import cocotb

from agent import ALUAgent
from scoreboard import ALUScoreboard
from coverage import ALUCoverage


class ALUEnv:

    def __init__(self, dut):

        self.agent = ALUAgent(dut)
        self.sb = ALUScoreboard()
        self.cov = ALUCoverage()

    def connect(self):

        self.agent.connect()

        self.sb.mon2sb = self.agent.m2sb
        self.cov.mon2cov = self.agent.m2cov

    async def run(self):

        cocotb.start_soon(self.agent.run())
        cocotb.start_soon(self.sb.run())
        cocotb.start_soon(self.cov.run())