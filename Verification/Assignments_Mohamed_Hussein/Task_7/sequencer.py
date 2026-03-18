import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

class Sequencer(uvm_sequencer):
    def __init__(self, name, parent):
        super().__init__(name, parent)