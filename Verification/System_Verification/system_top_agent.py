"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: system_top_agent.py

    Description:
        Top-level UVM agent for the full SYSTEM_TOP DUT.

        Active mode  (is_active = UVM_ACTIVE, default)
            Instantiates: sequencer + driver + monitor.
            The driver stimulates the DUT; the monitor observes all
            boundary signals.  Sequences are run through the sequencer.

        Passive mode  (is_active = UVM_PASSIVE)
            Instantiates: monitor only.
            Useful for black-box observation when another agent (or a
            physical interface) drives the DUT independently.

        The agent always exposes agt_ap (an analysis port that
        re-broadcasts whatever the monitor sees) so that a scoreboard
        or subscriber higher up the hierarchy can connect without
        needing to know whether the agent is active or passive.

        Configuration (set in the testbench env or test):
            ConfigDB().set(None, "uvm_test_top.*", "is_active", UVM_ACTIVE)
            ConfigDB().set(None, "uvm_test_top.*", "DUT",       dut)
"""

import pyuvm
from pyuvm import *

from system_top_sequencer import system_top_sequencer
from system_top_driver    import system_top_driver
from system_top_monitor   import system_top_monitor


class system_top_agent(uvm_agent):

    def build_phase(self):
        # ── Always-present analysis port ─────────────────────────────
        self.agt_ap = uvm_analysis_port("agt_ap", self)
        # ── Always-present monitor ───────────────────────────────────
        self.mon = system_top_monitor.create("mon", self)
        # ── Active-mode components (sequencer + driver) ──────────────
        self.sqr = system_top_sequencer.create("sqr", self)
        self.drv = system_top_driver.create("drv",    self)

    def connect_phase(self):
        # ── Always: re-broadcast monitor traffic up to the agent port ─
        self.mon.mon_ap.connect(self.agt_ap)
        # ── Active only: wire sequencer→driver ───────────────────────
        self.drv.seq_item_port.connect(self.sqr.seq_item_export)
