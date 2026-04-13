import random
from pyuvm import *


class ALU_seq_item(uvm_sequence_item):

    def __init__(self, name="ALU_seq_item"):
        super().__init__(name)

        self.a = 0
        self.b = 0
        self.op = 0
        self.out = 0
        self.c = 0

    # -------------------------------------------------
    # BASIC RANDOMIZATION
    # -------------------------------------------------
    def randomize(self):
        self.a = random.randint(0, 15)
        self.b = random.randint(0, 15)
        self.op = random.randint(0, 3)
        return self

    # -------------------------------------------------
    # CONSTRAINED RANDOMIZATION
    # -------------------------------------------------
    def randomize_with(self, constraint_fn):
        """
        Keep randomizing until constraint is satisfied
        """
        for _ in range(1000):
            self.randomize()
            if constraint_fn(self):
                return self

        raise ValueError("Constraint not satisfiable")

    # -------------------------------------------------
    # STRING PRINT
    # -------------------------------------------------
    def __str__(self):
        op_map = {0: "ADD", 1: "XOR", 2: "AND", 3: "OR"}
        op_str = op_map.get(self.op, "UNKNOWN")
        return f"{op_str} a={self.a} b={self.b} | out={self.out} c={self.c}"
