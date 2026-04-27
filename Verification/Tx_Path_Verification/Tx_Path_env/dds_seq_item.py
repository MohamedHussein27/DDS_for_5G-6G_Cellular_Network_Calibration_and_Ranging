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
        # Output variables 
        self.final_amplitude = 0

    def randomize(self):
        """Custom randomization mimicking SV constraints."""
        # 32-bit Tuning Word Start (Full range)
        self.FTW_start = random.randint(0, (1 << 32) - 1)
        
        # 13-bit Cycles (Ensure it is at least 1 to prevent immediate resets)
        self.cycles = random.randint(1, (1 << 13) - 1)
        
        # 32-bit Step. Constrained to 20 bits to prevent aliasing 
        # (stepping past the Nyquist limit in a single clock cycle).
        self.FTW_step = random.randint(0, (1 << 32) - 1)
        # Reset can be either active (0) or inactive (1)
        self.rst_n = random.choices([0, 1], weights=[10, 90])[0]  # Bias towards non-reset states for more meaningful tests
        
    
    # STRING FORMATTING
    def convert2string(self):
        """Full transaction printout (Inputs -> Outputs)"""
        return (f"{self.get_name()} "
                f"IN[ rst_n={self.rst_n}, FTW_start={hex(self.FTW_start)}, "
                f"FTW_step={hex(self.FTW_step)}, cycles={self.cycles} ] -> "
                f"OUT[ final_amplitude={hex(self.final_amplitude)} ]")

    def convert2string_stimulus(self):
        """Input-only printout for the driver/sequencer logging"""
        return (f"{self.get_name()} "
                f"STIMULUS: rst_n={self.rst_n}, FTW_start={hex(self.FTW_start)}, "
                f"FTW_step={hex(self.FTW_step)}, cycles={self.cycles}")       