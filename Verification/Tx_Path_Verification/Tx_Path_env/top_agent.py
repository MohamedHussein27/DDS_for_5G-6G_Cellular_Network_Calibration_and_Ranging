"""
"""

import pyuvm
from pyuvm import *
from top_sequencer import top_sequencer
from top_driver import top_driver
from top_monitor import top_monitor

class top_agent(uvm_agent):
    def build_phase(self):
        # 1. The Monitor is ALWAYS built, because we always want to observe
        self.mon = top_monitor.create("mon", self)
        self.agent_port = uvm_analysis_port("agent_port", self)
        self.is_active = ConfigDB().get(self,"","is_active")
        # 2. The Driver and Sequencer are ONLY built if the agent is ACTIVE
        if self.is_active == UVM_ACTIVE:
            self.sqr = top_sequencer.create("sqr", self)
            self.drv = top_driver.create("drv", self)

    def connect_phase(self):
        # 1. Broadcast the monitor's observed traffic up to the agent level
        self.mon.mon_port.connect(self.agent_port)
        
        # 2. Connect Sequencer to Driver ONLY if they were built
        if self.is_active == UVM_ACTIVE:
            self.drv.seq_item_port.connect(self.sqr.seq_item_export)