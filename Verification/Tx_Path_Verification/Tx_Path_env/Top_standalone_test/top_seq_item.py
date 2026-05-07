import random
import math 
from pyuvm import uvm_sequence_item

class top_item(uvm_sequence_item):
    def __init__(self, name="top_item"):
        super().__init__(name)
        
        # System Constants
        self.Fs = 491.52e6
        self.M  = 32  
        
        # High-Level Physical Parameters
        self.f0 = 0.0
        self.B  = 0.0

        # Hardware Inputs
        self.FTW_start = 0
        self.dds_enable = 0
        self.cycles = 4096
        self.FTW_step = 0
        self.rst_n = 1
        self.enable = 0
        
        # Output variables 
        self.final_amplitude = 0
        self.valid_out = 0
        
        # Streaming Metadata
        self.sample_index = 0 

    def calculate_hw_registers(self):
        """Calculates the generic equations based on the current f0, B, and cycles."""
        # FTW_start = (f0 * 2^M) / Fs
        self.FTW_start = int(math.floor((self.f0 * (1 << self.M)) / self.Fs))

        # FTW_step = (B * 2^M) / (Fs * Ns)
        self.FTW_step = int(math.floor((self.B * (1 << self.M)) / (self.Fs * self.cycles)))

        # Mask to ensure they strictly fit within M bits
        bit_mask = (1 << self.M) - 1
        self.FTW_start &= bit_mask
        self.FTW_step  &= bit_mask

    def set_fixed_chirp(self, f0, B, cycles=4096):
        """Helper function to cleanly force a specific chirp without randomizing."""
        self.f0 = f0
        self.B = B
        self.cycles = cycles
        self.rst_n = 1
        self.enable = 1
        # Trigger the math update
        self.calculate_hw_registers()

    def randomize(self):
        """Custom randomization for general testcases."""
        nyquist = self.Fs / 2
        self.f0 = random.uniform(0, nyquist / 2)             
        self.B  = random.uniform(1e6, nyquist - self.f0)     

        self.rst_n = random.choices([0, 1], weights=[5, 95])[0]
        self.enable = 1 

        # Trigger the math update using the newly randomized values
        self.calculate_hw_registers()

    def convert2string(self):
        return (f"{self.get_name()} "
                f"IN[ rst_n={self.rst_n}, f0={self.f0/1e6:.2f}MHz, B={self.B/1e6:.2f}MHz, "
                f"FTW_start={hex(self.FTW_start)}, FTW_step={hex(self.FTW_step)}, cycles={self.cycles} ] -> "
                f"OUT[ final_amplitude={hex(self.final_amplitude)} ]")

    def convert2string_stimulus(self):
        return (f"{self.get_name()} "
                f"STIMULUS: rst_n={self.rst_n}, f0={self.f0/1e6:.2f}MHz, B={self.B/1e6:.2f}MHz, "
                f"FTW_start={hex(self.FTW_start)}, FTW_step={hex(self.FTW_step)}, cycles={self.cycles}")