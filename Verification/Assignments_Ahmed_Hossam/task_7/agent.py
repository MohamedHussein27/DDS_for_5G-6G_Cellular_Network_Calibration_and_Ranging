import pyuvm
from pyuvm import *
from sequencer import sequencer
from driver import driver
from monitor import monitor

class agent(uvm_agent):
    def build_phase(self):
        
        self.sqr = sequencer("sequencer", self)
        self.drv = driver("driver", self)
        self.mon = monitor("monitor", self)
        
       
        self.agent_port = uvm_analysis_port("agent_port", self)

    def connect_phase(self):
        
        self.mon.mon_port.connect(self.agent_port)
        
        
        self.drv.seq_item_port.connect(self.sqr.seq_item_export)