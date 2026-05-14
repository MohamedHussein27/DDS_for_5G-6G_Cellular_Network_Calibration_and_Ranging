"""
===============================================================================
Sponsor      : Analog Devices, Inc. (ADI)
Institution  : Faculty of Engineering, Ain Shams University
Project      : DDS for 5G/6G Cellular Network Calibration and Ranging

Module       : agent.py

Description  :
    This module implements the DDS UVM Agent.

    The agent acts as a container for:
        - Sequencer
        - Driver
        - Monitor

    Responsibilities:
        1. Build active/passive components based on configuration.
        2. Connect driver and sequencer interfaces.
        3. Forward monitored transactions through an analysis port.

    Active Agent:
        Includes Driver + Sequencer + Monitor

    Passive Agent:
        Includes Monitor only

===============================================================================
"""

import pyuvm
from pyuvm import *

from dds_sequencer import dds_sequencer
from dds_driver import dds_driver
from dds_monitor import dds_monitor


class dds_agent(uvm_agent):

    def build_phase(self):
        """
        Build agent components depending on active/passive mode.
        """

        # ---------------------------------------------------------------------
        # Monitor is always created to observe DUT activity
        # ---------------------------------------------------------------------
        self.mon = dds_monitor.create("mon", self)

        # Agent-level analysis port used to broadcast monitored transactions
        self.agt_ap = uvm_analysis_port("agt_ap", self)

        # Retrieve active/passive configuration from ConfigDB
        self.is_active = ConfigDB().get(self, "", "is_active")

        # ---------------------------------------------------------------------
        # Build Driver and Sequencer only for ACTIVE agents
        # ---------------------------------------------------------------------
        if self.get_is_active() == uvm_active_passive_enum.UVM_ACTIVE:

            # Sequencer generates transaction flow
            self.sqr = dds_sequencer.create("sqr", self)

            # Driver converts transactions into DUT pin activity
            self.drv = dds_driver.create("drv", self)

    def connect_phase(self):
        """
        Connect agent internal communication ports.
        """

        # Forward monitor transactions to upper verification layers
        self.mon.mon_ap.connect(self.agt_ap)
        
        # ---------------------------------------------------------------------
        # Connect Sequencer <-> Driver only in ACTIVE mode
        # ---------------------------------------------------------------------
        if self.get_is_active() == uvm_active_passive_enum.UVM_ACTIVE:

            # Driver pulls transactions from sequencer
            self.drv.seq_item_port.connect(self.sqr.seq_item_export)