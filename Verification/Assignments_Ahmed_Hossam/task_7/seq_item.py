import pyuvm
from pyuvm import *
import random
class seq_item (uvm_sequence_item):
    def __init__(self, name="seq_item"):
        super().__init__(name)
        self.a=0
        self.b=0
        self.op=0
        self.c=0
        self.out=0
        self.reset = 0
     
     
    def randomize(self):
        self.a = random.randint(0, 15)
        self.b = random.randint(0, 15)
        self.op = random.randint(0, 3)    