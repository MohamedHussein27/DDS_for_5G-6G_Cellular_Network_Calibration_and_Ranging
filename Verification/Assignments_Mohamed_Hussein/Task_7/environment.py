import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

from agent import *
from scoreboard import *
from subscriber import *

class Environment(uvm_env):
    def __init__(self, name, parent):
        super().__init__(name, parent)


    def build_phase(self):
        self.sb = Scoreboard.create("sb", self)
        self.sub = Subscriber.create("sub", self)
        self.agt = Agent.create("agt", self)


    def connect_phase(self):
        self.agt.agt_ap.connect(self.sb.sb_export)
        self.agt.agt_ap.connect(self.sub.analysis_export)