from pyuvm import *
from ALU_driver import ALU_driver
from ALU_monitor import ALU_monitor
from ALU_sequencer import ALU_sequencer


class ALU_agent(uvm_agent):
    def build_phase(self):
        super().build_phase()
        self.agent_ap = uvm_analysis_port("agent_ap", self)

        self.sr_ag = ALU_sequencer("sr_ag", self)
        self.dr_ag = ALU_driver("dr_ag", self)
        self.mon_ag = ALU_monitor("mon_ag", self)

    def connect_phase(self):
        super().connect_phase()
        self.dr_ag.seq_item_port.connect(self.sr_ag.seq_item_export)
        self.mon_ag.monitor_ap.connect(self.agent_ap)
