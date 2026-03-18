import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

# multiple inheritance
class Sequence_Item(uvm_sequence_item, Randomized):

    def __init__(self, name="Sequence_Item"):
        # 3. You must initialize both parent classes
        uvm_sequence_item.__init__(self, name)
        Randomized.__init__(self)

        self.rst_n = 1
        self.a = 0
        self.b = 0
        self.op = 0

        # Output variables (Not randomized)
        self.out = 0
        self.c = 0

        # =========================================================
        # CRV SETUP: Registering variables with the solver
        # =========================================================
        # Syntax: self.add_rand("variable_name", [list_of_valid_values])
        
        self.add_rand("rst_n", [0, 1])
        self.add_rand("a", list(range(16)))  # 4-bit (0 to 15)
        self.add_rand("b", list(range(16)))  # 4-bit (0 to 15)
        self.add_rand("op", list(range(4)))  # 2-bit (0 to 3)

    # =========================================================
    # STRING FORMATTING
    # =========================================================
    def convert2string(self):
        # The :b format specifier converts the integers to binary automatically
        return (f"{self.get_name()} rst_n = 0b{self.rst_n:b}, "
                f"A = 0b{self.a:b}, B = 0b{self.b:b}, OP = 0b{self.op:b}, "
                f"OUT = 0b{self.out:b}, Carry = 0b{self.c:b}")

    def convert2string_stimulus(self):
        return (f"{self.get_name()} rst_n = 0b{self.rst_n:b}, "
                f"A = 0b{self.a:b}, B = 0b{self.b:b}, OP = 0b{self.op:b}")