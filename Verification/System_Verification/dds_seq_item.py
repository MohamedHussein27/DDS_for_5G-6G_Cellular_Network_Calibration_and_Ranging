import random
from pyuvm import uvm_sequence_item

class dds_seq_item(uvm_sequence_item):
    def __init__(self, name="dds_seq_item"):
        super().__init__(name)
        # Hardware Inputs
        self.FTW_start = 0
        self.cycles = 0
        self.FTW_step = 0
        self.rst_n = 1
        self.enable = 0
        
        # Output variables 
        self.final_amplitude = 0
        self.valid_out = 0
        
        # Streaming Metadata
        self.sample_index = 0 # NEW: Tracks the current clock cycle of the stream

    def randomize(self):
        """Custom randomization mimicking SV constraints."""
        self.FTW_start = random.randint(0, (1 << 32) - 1)
        self.FTW_step = random.randint(0, (1 << 32) - 1)
        self.cycles = random.randint(0,(1 << 13) - 1)  
        self.rst_n = random.choices([0, 1], weights=[5, 95])[0]  
    
    def convert2string(self):
        return (f"{self.get_name()} "
                f"IN[ rst_n={self.rst_n}, FTW_start={hex(self.FTW_start)}, "
                f"FTW_step={hex(self.FTW_step)}, cycles={self.cycles} ] -> "
                f"OUT[ final_amplitude={hex(self.final_amplitude)} ]")

    def convert2string_stimulus(self):
        return (f"{self.get_name()} "
                f"STIMULUS: rst_n={self.rst_n}, FTW_start={hex(self.FTW_start)}, "
                f"FTW_step={hex(self.FTW_step)}, cycles={self.cycles}")