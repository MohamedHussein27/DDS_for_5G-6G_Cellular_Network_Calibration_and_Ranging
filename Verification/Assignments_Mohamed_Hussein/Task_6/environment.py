import cocotb
from cocotb.triggers import *
from agent import Agent
from scoreboard import Scoreboard
from subscriber import Subscriber

class Environment:

    def __init__(self):
        
        self.join_any = Event()

        self.agent = Agent()
        self.s = Scoreboard()
        self.sub = Subscriber()

    def connect(self):

        self.agent.connect()

        self.s.mon2sb = self.agent.mon2sb
        self.sub.mon2cov = self.agent.mon2cov

        # finishing simulation
        self.agent.gen.finish = self.join_any

    async def run_environment(self, dut, test_type):

        cocotb.start_soon(self.agent.run(dut, test_type))
        cocotb.start_soon(self.s.run())
        cocotb.start_soon(self.sub.run())