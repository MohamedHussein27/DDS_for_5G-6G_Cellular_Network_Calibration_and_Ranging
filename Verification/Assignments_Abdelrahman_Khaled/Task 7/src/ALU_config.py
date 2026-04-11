from pyuvm import *
import cocotb
from cocotb.triggers import RisingEdge, FallingEdge

# ALU_config.sv converted


class ALU_config(uvm_object):
    def __init__(self, name="ALU_config"):
        super().__init__(name)
        self.vif = None  # This will hold the cocotb DUT handle
