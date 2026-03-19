from pyuvm import *
import random
from alu_pkg import opcode_t

class alu_sequence_item(uvm_sequence_item):
    def __init__(self, name="alu_sequence_item"):
        super().__init__(name)
        self.a = 0
        self.b = 0
        self.op = opcode_t.ADD_op
        self.rst_n = 1
        self.result = 0

    def randomize(self):
        self.a = random.randint(0, 15)
        self.b = random.randint(0, 15)
        self.op = random.choice(list(opcode_t))
        self.rst_n = random.choices([0, 1], weights=[1, 100], k=1)[0]
        return True

    def convert2string(self):
        # Prevents crashes if a sequence passes a raw integer
        try:
            op_name = opcode_t(self.op).name
        except ValueError:
            op_name = "UNKNOWN"
            
        return f"A={self.a:2} B={self.b:2} OP={op_name:6} RST_N={self.rst_n} RES={self.result:2}"