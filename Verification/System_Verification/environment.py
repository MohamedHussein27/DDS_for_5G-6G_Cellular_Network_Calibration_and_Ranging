"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
    
    Module: environment.py

    Description:
        The Environment class acts as the top-level container for the TX wrapper 
        testbench. Following the verification architecture, it instantiates the 
        Top Active Agent to drive the DUT, alongside three Passive Agents (DDS, 
        FFT, IFFT) to monitor internal DSP boundaries.
"""

import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import logging

# 1. Top Component Imports
from rx_top_agent import rx_top_agent
from rx_top_scoreboard import *
from rx_top_subscriber import *


from tx_top_agent import tx_top_agent
from tx_top_scoreboard import *
from tx_top_subscriber import *


from system_top_agent import system_top_agent
from system_top_scoreboard import *
from system_top_subscriber import *


class Environment(uvm_env):
    def __init__(self, name, parent):
        super().__init__(name, parent)

    def build_phase(self):
        super().build_phase()

        # 1. Read the configuration (Default to "TOP" if not set)
        if ConfigDB().exists(self, "", "VERIF_MODE"):
            self.mode = ConfigDB().get(self, "", "VERIF_MODE")
        else:
            self.mode = "SYSTEM"

        self.logger.info(f"Environment building in {self.mode} mode.")

        # 2. Conditionally Build Components (switch-case style)

        # ---------------------------------------------------------
        # TX PATH 
        # ---------------------------------------------------------
        if self.mode == "TX":
            self.tx_agt = tx_top_agent.create("tx_agt", self)
            self.tx_sb  = tx_top_scoreboard_fp.create("tx_sb", self)
            self.tx_sub = tx_top_subscriber.create("tx_sub", self)
            
        
        # ---------------------------------------------------------
        # FFT PATH 
        # ---------------------------------------------------------
        elif self.mode == "RX":
            self.rx_agt = rx_top_agent.create("rx_agt", self)
            self.rx_sb  = rx_top_scoreboard.create("rx_sb", self)
            self.rx_sub = rx_top_subscriber.create("rx_sub", self)
            
        # ---------------------------------------------------------
        # SYSTEM PATH 
        # ---------------------------------------------------------
        elif self.mode == "SYSTEM":
            # Build EVERYTHING for system-level regression
            
            # Top Active Path
            self.system_agt  = system_top_agent.create("system_agt", self)
            self.system_sb   = system_top_scoreboard.create("system_sb", self)
            self.system_sub  = system_top_subscriber.create("system_sub", self)
            
            # IFFT Passive Path
            self.tx_agt = tx_top_agent.create("tx_agt", self)
            self.tx_sb  = tx_top_scoreboard_fp.create("tx_sb", self)
            self.tx_sub = tx_top_subscriber.create("tx_sub", self)
            
            
            # FFT Passive Path
            self.rx_agt = rx_top_agent.create("rx_agt", self)
            self.rx_sb  = rx_top_scoreboard.create("rx_sb", self)
            self.rx_sub = rx_top_subscriber.create("rx_sub", self)
            
            
        else:
            self.logger.fatal(f"Unknown VERIF_MODE: {self.mode}")


    def connect_phase(self):
        super().connect_phase()

        
        if self.mode == "SYSTEM":
            # Connect Top Agent to its Scoreboard and Subscriber
            self.system_agt.agt_ap.connect(self.system_sb.sb_export)
            self.system_agt.agt_ap.connect(self.system_sub.analysis_export)

            
            # Connect FFT Agent to its Scoreboard and Subscriber
            self.tx_agt.agt_ap.connect(self.tx_sb.sb_export)
            self.tx_agt.agt_ap.connect(self.tx_sub.analysis_export)

            # Connect IFFT Agent to its Scoreboard and Subscriber
            self.rx_agt.agt_ap.connect(self.rx_sb.sb_export)
            self.rx_agt.agt_ap.connect(self.rx_sub.analysis_export)
        
      

        elif self.mode == "TX":
            # Connect FFT Agent to its Scoreboard and Subscriber
            self.tx_agt.agt_ap.connect(self.tx_sb.sb_export)
            self.tx_agt.agt_ap.connect(self.tx_sub.analysis_export)
        
        elif self.mode == "RX":
            # Connect IFFT Agent to its Scoreboard and Subscriber
            self.rx_agt.agt_ap.connect(self.rx_sb.sb_export)
            self.rx_agt.agt_ap.connect(self.rx_sub.analysis_export)