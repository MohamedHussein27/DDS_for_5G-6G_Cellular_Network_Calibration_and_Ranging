import pyuvm
from pyuvm import *
from seq_item import *
from cocotb.triggers import Timer
from cocotb.clock import Clock
from driver import *
from subscriber import subscriber
from scoreboard import scoreboard
from agent import agent
class env(uvm_env):
    def build_phase (self):
        self.agent_ = agent.create("agent", self)
        self.scoreboard = scoreboard.create("scoreboard", self)
        self.subscriber = subscriber.create("subscriber", self)

    def connect_phase(self): 
        self.agent_.agent_port.connect(self.scoreboard.sc_export)
        self.agent_.agent_port.connect(self.subscriber.sub_export)   