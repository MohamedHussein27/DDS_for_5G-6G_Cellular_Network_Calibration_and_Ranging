"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_agent.py

    Description:
        The IFFT agent bundles the sequencer, driver, and monitor into a single 
        hierarchical block. Based on the UVM configuration database, it acts as 
        an Active Agent (injecting stimulus) for block-level testing, or a Passive 
        Agent (monitoring only) for top-level TX path integration.
"""

from pyuvm import *
from ifft_driver import ifft_driver
from ifft_monitor import ifft_monitor


class ifft_agent(uvm_agent):
    def build_phase(self):
        # 1. The Monitor is ALWAYS built to observe traffic
        self.mon = ifft_monitor.create("mon", self)
        self.agt_ap = uvm_analysis_port("agt_ap", self)

        # getting the active/passive configuration from the Environment/Base Test
        self.is_active = ConfigDB().get(self,"","is_active")

        # 2. Driver and Sequencer are ONLY built if the agent is ACTIVE
        if self.is_active == uvm_active_passive_enum.UVM_ACTIVE:
            self.sqr = uvm_sequencer.create("sqr", self)
            self.drv = ifft_driver.create("drv", self)

    def connect_phase(self):
        # 1. Broadcast the monitor's observed traffic up to the agent boundary
        self.mon.mon_port.connect(self.agt_ap)
        
        # 2. Connect the Sequencer to the Driver ONLY if they were built
        if self.is_active == uvm_active_passive_enum.UVM_ACTIVE:
            self.drv.seq_item_port.connect(self.sqr.seq_item_export)