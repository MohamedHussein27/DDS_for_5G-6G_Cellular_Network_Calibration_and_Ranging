import pyuvm
from pyuvm import *
from seq_item import *
from driver import *
from sequencer import sequencer
from monitor import monitor
from cocotb.triggers import Timer
from cocotb.clock import Clock
class agent(uvm_agent):
    def build_phase(self):
        self.sequencer = sequencer.create("sequencer", self)
        self.driver = driver.create("driver", self)
        self.monitor = monitor.create("monitor", self)
        self.agent_port = uvm_analysis_port("agent_port",self)
        self.dut = ConfigDB().get(self, "", "DUT")
        
    def connect_phase(self):
        self.monitor.mon_port.connect(self.agent_port)
        self.driver.seq_item_port.connect(self.sequencer.seq_item_export)   
        self.driver.dut = self.dut
        self.monitor.dut = self.dut