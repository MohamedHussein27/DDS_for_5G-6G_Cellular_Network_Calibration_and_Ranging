from pyuvm import *
from fft_seq_item import fft_item

class FftSequencer(uvm_sequencer):
    def __init__(self, name="fft_sequencer", parent=None):
        super().__init__(name, parent)