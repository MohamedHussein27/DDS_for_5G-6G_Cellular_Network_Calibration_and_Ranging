import pyuvm
from pyuvm import *
from alu_agent import alu_agent
from alu_scoreboard import alu_scoreboard
from alu_coverage import alu_coverage

class alu_env(uvm_env):
    def build_phase(self):
        super().build_phase()
        self.agt = alu_agent("agt", self)
        self.sb = alu_scoreboard("sb", self)
        self.cov = alu_coverage("cov", self)

    def connect_phase(self):
        super().connect_phase()
        # Ports use .connect() to attach to the Exports we assigned in the children
        self.agt.agt_ap.connect(self.sb.sb_export)
        self.agt.agt_ap.connect(self.cov.analysis_export)