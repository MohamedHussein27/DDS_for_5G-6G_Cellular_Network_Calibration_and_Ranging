"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: rx_agent.py  (system-level wrapper)

    Description:
        Active/passive RX sub-agent for use inside the system-level
        UVM environment (system_top_env).

        It wraps the existing RX components (top_sequencer, rx_driver,
        rx_monitor) that were originally developed for the RX_TOP
        standalone testbench.  By re-hosting them under a configurable
        agent, the system environment can optionally enable an
        independent RX stimulus path (e.g. for injection tests) without
        changing the top-level test hierarchy.

        Active mode  (UVM_ACTIVE, default)
            Instantiates: top_sequencer + rx_driver + rx_monitor.
            Sequences can independently drive the RX datapath through
            env.rx_agent.sqr.seq_item_export.

        Passive mode  (UVM_PASSIVE)
            Instantiates: rx_monitor only.
            The loopback wires from TX feed RX automatically inside
            SYSTEM_TOP, so in most system tests the RX agent is passive
            and only observes.

        The agent always exposes agt_ap so the system scoreboard and
        subscriber can tap RX transactions independently.

        NOTE: In SYSTEM_TOP, RX inputs (rx_valid_in, rx_in_re, rx_in_im,
        ref_wr_en, ref_wr_re, ref_wr_im) are driven internally by the
        TX output wires.  When the RX agent is in ACTIVE mode it can
        override those wires — use with care to avoid bus contention.
"""

import pyuvm
from pyuvm import *

# Re-use the RX standalone components directly
from rx_top_sequencer import top_sequencer   # shared sequencer class
from rx_top_driver  import rx_driver
from rx_top_monitor import rx_monitor


class rx_top_agent(uvm_agent):

    def build_phase(self):
        # ── Always-present analysis port ─────────────────────────────
        self.agt_ap = uvm_analysis_port("agt_ap", self)

        # ── Always-present monitor ───────────────────────────────────
        self.mon = rx_monitor.create("mon", self)

        # ── Active-mode only ─────────────────────────────────────────
        try:
            self.is_active = ConfigDB().get(self, "", "is_active")
        except UVMConfigItemNotFound:
            self.is_active = uvm_active_passive_enum.UVM_ACTIVE  # Default to active for block-level

        # 3. Driver and Sequencer are ONLY built if the agent is ACTIVE
        if self.is_active == uvm_active_passive_enum.UVM_ACTIVE:
            self.sqr = top_sequencer("sqr", self)
            self.drv = rx_driver("drv", self)

    def connect_phase(self):
        # Always re-broadcast monitor observations to the agent port
        self.mon.mon_ap.connect(self.agt_ap)

        # Active only: connect sequencer → driver
        if self.is_active == uvm_active_passive_enum.UVM_ACTIVE:
            self.drv.seq_item_port.connect(self.sqr.seq_item_export)