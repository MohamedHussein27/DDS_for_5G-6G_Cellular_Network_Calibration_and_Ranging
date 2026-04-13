from pyuvm import *
from ALU_agent import ALU_agent
from ALU_scoreboard import ALU_scoreboard


class ALU_env(uvm_env):
    def build_phase(self):
        super().build_phase()
        self.ag_env = ALU_agent("ag_env", self)
        self.sb_env = ALU_scoreboard("sb_env", self)

    def connect_phase(self):
        super().connect_phase()
        # Connect monitor to scoreboard
        self.ag_env.agent_ap.connect(self.sb_env.sb_export.analysis_export)
