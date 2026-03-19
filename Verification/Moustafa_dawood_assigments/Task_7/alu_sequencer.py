from pyuvm import *

class alu_sequencer(uvm_sequencer):
    def __init__(self, name, parent):
        super().__init__(name, parent)