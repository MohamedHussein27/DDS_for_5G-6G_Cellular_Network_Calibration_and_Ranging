"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: sequencer.py

    Description:
        This module defines the PUVM Sequencer for the DDS stimulus generation.
        It acts as the traffic controller (conveyor belt) between the high-level 
        test sequences and the pin-level Driver.

        Key responsibilities:
        1. Arbitration: Manages requests from multiple sequences (e.g., Frequency 
           Tuning Word updates, Phase offsets, and Reset sequences like Group F).
        2. Routing: Passes the randomized `seq_item` transaction packets downstream 
           to the Driver precisely when the hardware is ready to accept them.
"""

import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

class Sequencer(uvm_sequencer):
    pass