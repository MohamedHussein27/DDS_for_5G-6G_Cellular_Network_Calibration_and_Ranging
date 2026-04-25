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
from top_agent import top_agent
from top_scoreboard import top_scoreboard
from top_subscriber import top_subscriber

# 2. DDS Component Imports (Added missing imports)
from dds_agent import dds_agent
from dds_scoreboard import dds_scoreboard
from dds_subscriber import dds_subscriber

# 3. FFT Component Imports
from fft_agent import fft_agent
from fft_scoreboard import fft_scoreboard
from fft_subscriber import fft_subscriber

# 4. IFFT Component Imports
from ifft_agent import ifft_agent
from ifft_scoreboard import ifft_scoreboard
from ifft_subscriber import ifft_subscriber


class Environment(uvm_env):
    def __init__(self, name, parent):
        super().__init__(name, parent)

    def build_phase(self):
        super().build_phase()

        # 1. Read the configuration (Default to "TOP" if not set)
        if ConfigDB().exists(self, "", "VERIF_MODE"):
            self.mode = ConfigDB().get(self, "", "VERIF_MODE")
        else:
            self.mode = "TOP"

        self.logger.info(f"Environment building in {self.mode} mode.")

        # 2. Conditionally Build Components (switch-case style)

        # ---------------------------------------------------------
        # IFFT PATH 
        # ---------------------------------------------------------
        if self.mode == "IFFT":
            self.ifft_agt = ifft_agent.create("ifft_agt", self)
            self.ifft_sb  = ifft_scoreboard.create("ifft_sb", self)
            self.ifft_sub = ifft_subscriber.create("ifft_sub", self)
            
        # ---------------------------------------------------------
        # DDS PATH 
        # ---------------------------------------------------------
        elif self.mode == "DDS":
            self.dds_agt = dds_agent.create("dds_agt", self)
            self.dds_sb  = dds_scoreboard.create("dds_sb", self)
            self.dds_sub = dds_subscriber.create("dds_sub", self)
        
        # ---------------------------------------------------------
        # FFT PATH 
        # ---------------------------------------------------------
        elif self.mode == "FFT":
            self.fft_agt = fft_agent.create("fft_agt", self)
            self.fft_sb  = fft_scoreboard.create("fft_sb", self)
            self.fft_sub = fft_subscriber.create("fft_sub", self)
            
        # ---------------------------------------------------------
        # TOP PATH 
        # ---------------------------------------------------------
        elif self.mode == "TOP":
            # Build EVERYTHING for system-level regression
            
            # Top Active Path
            self.top_agt  = top_agent.create("top_agt", self)
            self.top_sb   = top_scoreboard.create("top_sb", self)
            self.top_sub  = top_subscriber.create("top_sub", self)
            
            # IFFT Passive Path
            self.ifft_agt = ifft_agent.create("ifft_agt", self)
            self.ifft_sb  = ifft_scoreboard.create("ifft_sb", self)
            self.ifft_sub = ifft_subscriber.create("ifft_sub", self)
            
            # DDS Passive Path
            self.dds_agt  = dds_agent.create("dds_agt", self)
            self.dds_sb   = dds_scoreboard.create("dds_sb", self)
            self.dds_sub  = dds_subscriber.create("dds_sub", self)
            
            # FFT Passive Path
            self.fft_agt  = fft_agent.create("fft_agt", self)
            self.fft_sb   = fft_scoreboard.create("fft_sb", self)
            self.fft_sub  = fft_subscriber.create("fft_sub", self)
            
        else:
            self.logger.fatal(f"Unknown VERIF_MODE: {self.mode}")


    def connect_phase(self):
        super().connect_phase()

        
        if self.mode == "TOP":
            # Connect Top Agent to its Scoreboard and Subscriber
            self.top_agt.agt_ap.connect(self.top_sb.sb_export)
            self.top_agt.agt_ap.connect(self.top_sub.analysis_export)

            # Connect DDS Agent to its Scoreboard and Subscriber
            self.dds_agt.agt_ap.connect(self.dds_sb.sb_export)
            self.dds_agt.agt_ap.connect(self.dds_sub.analysis_export)

            # Connect FFT Agent to its Scoreboard and Subscriber
            self.fft_agt.agt_ap.connect(self.fft_sb.sb_export)
            self.fft_agt.agt_ap.connect(self.fft_sub.analysis_export)

            # Connect IFFT Agent to its Scoreboard and Subscriber
            self.ifft_agt.agt_ap.connect(self.ifft_sb.sb_export)
            self.ifft_agt.agt_ap.connect(self.ifft_sub.analysis_export)
        
        elif self.mode == "DDS":
            # Connect DDS Agent to its Scoreboard and Subscriber
            self.dds_agt.agt_ap.connect(self.dds_sb.sb_export)
            self.dds_agt.agt_ap.connect(self.dds_sub.analysis_export)

        elif self.mode == "FFT":
            # Connect FFT Agent to its Scoreboard and Subscriber
            self.fft_agt.agt_ap.connect(self.fft_sb.sb_export)
            self.fft_agt.agt_ap.connect(self.fft_sub.analysis_export)
        
        elif self.mode == "IFFT":
            # Connect IFFT Agent to its Scoreboard and Subscriber
            self.ifft_agt.agt_ap.connect(self.ifft_sb.sb_export)
            self.ifft_agt.agt_ap.connect(self.ifft_sub.analysis_export)