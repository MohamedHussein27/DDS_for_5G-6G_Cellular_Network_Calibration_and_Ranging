import pyuvm
from pyuvm import *
from alu_sequencer import alu_sequencer
from alu_driver import alu_driver
from alu_monitor import alu_monitor

class alu_agent(uvm_agent):
    def build_phase(self):
        super().build_phase()
        self.sqr = alu_sequencer("sqr", self)
        self.driver = alu_driver("driver", self)
        self.monitor = alu_monitor("monitor", self)
        self.agt_ap = uvm_analysis_port("agt_ap", self)
        
        # Standard UVM retrieval matching the wildcard "*" from the test
        self.dut = ConfigDB().get(self, "", "DUT")

    def connect_phase(self):
        super().connect_phase()
        self.driver.seq_item_port.connect(self.sqr.seq_item_export)
        self.monitor.mon_ap.connect(self.agt_ap)
        self.driver.dut = self.dut
        self.monitor.dut = self.dut