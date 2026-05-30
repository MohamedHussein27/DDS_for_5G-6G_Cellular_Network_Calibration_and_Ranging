"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: system_top_sequencer.py

    Description:
        Sequencer for the SYSTEM_TOP agent.
        Passes system_top_item transactions to the system driver.
        Mirrors the pattern of top_sequencer and rx_top_sequencer.
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
import logging


class system_top_sequencer(uvm_sequencer):
    pass
