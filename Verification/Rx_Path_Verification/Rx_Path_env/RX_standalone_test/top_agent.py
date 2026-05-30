"""
"""

import pyuvm
from pyuvm import *
from top_sequencer import top_sequencer
from top_driver import rx_driver
from top_monitor import rx_monitor

class top_agent(uvm_agent):
    def build_phase(self):       
        self.agt_ap = uvm_analysis_port("agt_ap", self)
        self.sqr = top_sequencer.create("sqr", self)
        self.drv = rx_driver.create("drv", self)
        self.mon = rx_monitor.create("mon", self)

    def connect_phase(self):
        # 1. Broadcast the monitor's observed traffic up to the agent level
        self.mon.mon_ap.connect(self.agt_ap)
        
        # 2. Connect Sequencer to Driver 
        self.drv.seq_item_port.connect(self.sqr.seq_item_export)