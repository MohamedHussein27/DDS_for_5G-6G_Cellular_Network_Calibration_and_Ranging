import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

from monitor import *
from driver import *
from sequencer import *


class Agent(uvm_agent):
    def __init__(self, name, parent):
        super().__init__(name, parent)


    def build_phase(self):
        # agent analysis port
        self.agt_ap = uvm_analysis_port.create("agt_ap", self)

        self.mon = Monitor.create("mon", self)
        self.drv = Driver.create("drv", self)
        self.sqr = Sequencer.create("sqr", self)

        


    def connect_phase(self):
        # connecting sequencer and driver
        self.drv.seq_item_port.connect(self.sqr.seq_item_export)

        # connecting agent port with monitor port
        self.mon.mon_ap.connect(self.agt_ap)
        