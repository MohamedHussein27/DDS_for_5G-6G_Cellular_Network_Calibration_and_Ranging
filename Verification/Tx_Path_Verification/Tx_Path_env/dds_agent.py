"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: agent.py

    Description:
        This module defines the UVM Agent for the DDS TX datapath.
        It serves as the top-level container and configuration hub for the 
        component triad: Sequencer, Driver, and Monitor. 
        
        Key responsibilities:
        1. Configuration: Evaluates the `is_active` flag to determine if the 
           environment requires active stimulus (Driver + Sequencer) or 
           passive observation (Monitor only).
        2. Instantiation: Builds the required sub-components during the build phase.
        3. Connection: Wires the Sequencer's export to the Driver's port, and 
           exposes the Monitor's analysis port to the broader environment/scoreboard.
"""

import pyuvm
from pyuvm import *
from sequencer import sequencer
from driver import driver
from monitor import monitor

class agent(uvm_agent):
    def build_phase(self):
        # 1. The Monitor is ALWAYS built, because we always want to observe
        self.mon = monitor.create("monitor", self)
        self.agent_port = uvm_analysis_port("agent_port", self)
        self.is_active = ConfigDB().get(self,"","is_active")
        # 2. The Driver and Sequencer are ONLY built if the agent is ACTIVE
        if self.is_active == UVM_ACTIVE:
            self.sqr = sequencer.create("sequencer", self)
            self.drv = driver.create("driver", self)

    def connect_phase(self):
        # 1. Broadcast the monitor's observed traffic up to the agent level
        self.mon.mon_port.connect(self.agent_port)
        
        # 2. Connect Sequencer to Driver ONLY if they were built
        if self.is_active == UVM_ACTIVE:
            self.drv.seq_item_port.connect(self.sqr.seq_item_export)