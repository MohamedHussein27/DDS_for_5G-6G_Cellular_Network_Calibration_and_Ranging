"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: tx_agent.py  (system-level wrapper)

    Description:
        Active/passive TX sub-agent for use inside the system-level
        UVM environment (system_top_env).

        It wraps the existing TX components (top_sequencer, top_driver,
        top_monitor) that were originally developed for the TX_TOP
        standalone testbench.  By re-hosting them inside a configurable
        agent, the system environment can turn TX stimulus on/off
        without restructuring the test hierarchy.

        Active mode  (UVM_ACTIVE, default)
            Instantiates: top_sequencer + top_driver + top_monitor.
            Sequences access the TX DUT through
            env.tx_agent.sqr.seq_item_export.

        Passive mode  (UVM_PASSIVE)
            Instantiates: top_monitor only.
            Useful for system-integration scenarios where TX is
            driven by another master.

        The agent always exposes agt_ap so the system scoreboard and
        subscriber can tap TX transactions.

        NOTE: The TX driver uses self.dut.u_tx for backdoor ROM writes
        ("u_tx" is the TX_TOP instance name inside SYSTEM_TOP).
        Make sure the testbench registers the DUT with key "DUT"
        pointing to the SYSTEM_TOP handle before run_phase begins.
"""

import pyuvm
from pyuvm import *

# Re-use the TX standalone components directly
from tx_top_sequencer import top_sequencer
from tx_top_driver    import top_driver
from tx_top_monitor   import top_monitor


class tx_top_agent(uvm_agent):

    def build_phase(self):
        # ── Always-present analysis port ─────────────────────────────
        self.agt_ap = uvm_analysis_port("agt_ap", self)

        # ── Always-present monitor ───────────────────────────────────
        self.mon = top_monitor.create("mon", self)

        # ── Active-mode only ─────────────────────────────────────────
        try:
            self.is_active = ConfigDB().get(self, "", "is_active")
        except UVMConfigItemNotFound:
            self.is_active = uvm_active_passive_enum.UVM_ACTIVE  # Default to active for block-level

        # 3. Driver and Sequencer are ONLY built if the agent is ACTIVE
        if self.is_active == uvm_active_passive_enum.UVM_ACTIVE:
            self.sqr = top_sequencer("sqr", self)
            self.drv = top_driver("drv", self)

    def connect_phase(self):
        # Always re-broadcast monitor observations to the agent port
        self.mon.mon_ap.connect(self.agt_ap)

        # Active only: connect sequencer → driver
        if self.is_active == uvm_active_passive_enum.UVM_ACTIVE:
            self.drv.seq_item_port.connect(self.sqr.seq_item_export)