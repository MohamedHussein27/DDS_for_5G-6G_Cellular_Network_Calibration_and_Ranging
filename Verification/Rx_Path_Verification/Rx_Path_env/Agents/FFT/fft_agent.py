"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: fft_agent.py

    Description:
        The FFT agent bundles the sequencer, driver, and monitor into a single 
        hierarchical block. Based on the UVM configuration database, it acts as 
        an Active Agent (injecting stimulus) for block-level testing, or a Passive 
        Agent (monitoring only) for top-level TX path integration.
"""

from pyuvm import *
from fft_driver import fft_driver
from fft_monitor import fft_monitor

class fft_agent(uvm_agent):
    def build_phase(self):
        super().build_phase()
        
        # 1. The Monitor is ALWAYS built to observe traffic
        self.mon = fft_monitor("mon", self)
        self.agt_ap = uvm_analysis_port("agt_ap", self)

        # 2. Get the active/passive configuration from the Environment/Base Test
        try:
            self.is_active = ConfigDB().get(self, "", "is_active")
        except UVMConfigItemNotFound:
            self.is_active = uvm_active_passive_enum.UVM_ACTIVE  # Default to active for block-level

        # 3. Driver and Sequencer are ONLY built if the agent is ACTIVE
        if self.is_active == uvm_active_passive_enum.UVM_ACTIVE:
            self.sqr = uvm_sequencer("sqr", self)
            self.drv = fft_driver("drv", self)

    def connect_phase(self):
        super().connect_phase()
        
        # 1. Broadcast the monitor's observed traffic up to the agent boundary
        self.mon.mon_port.connect(self.agt_ap)
        
        # 2. Connect the Sequencer to the Driver ONLY if they were built
        if self.is_active == uvm_active_passive_enum.UVM_ACTIVE:
            self.drv.seq_item_port.connect(self.sqr.seq_item_export)