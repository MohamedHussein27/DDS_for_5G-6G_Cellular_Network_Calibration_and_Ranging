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
        self.agent_ = agent("agent", self) 
        self.scoreboard = scoreboard("scoreboard", self)
        self.subscriber = subscriber("subscriber", self)

    def connect_phase(self):
        # 1. Route data to the Scoreboard
        self.agent_.agent_port.connect(self.scoreboard.sc_fifo.analysis_export)
        
        # 2. Route the exact same data to the Subscriber!
        self.agent_.agent_port.connect(self.subscriber.sub_fifo.analysis_export)