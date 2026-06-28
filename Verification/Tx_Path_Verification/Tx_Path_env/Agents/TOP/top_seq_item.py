import random
import math
from pyuvm import uvm_sequence_item

class top_item(uvm_sequence_item):
    def __init__(self, name="top_item"):
        super().__init__(name)
        
        # System Constants
        self.Fs = 491520000  # Exact integer
        self.M  = 32  
        
        # High-Level Physical Parameters
        self.f0 = 0
        self.B  = 0
        self.target_cycles = 4096

        # ==========================================
        # NEW: Physical Bus Inputs (To Driver)
        # ==========================================
        self.rst_n = 1
        self.addr = 0
        self.wr_en = 0
        self.wr_data = 0
        self.rd_en = 0
        
        # Internal Calculated Registers (No longer driven directly)
        self.FTW_start = 0
        self.FTW_step = 0
        
        # Output variables 
        self.rd_data = 0
        self.dds_ready_flag = 0
        self.final_amplitude = 0
        self.valid_out = 0
        self.sample_index = 0 

        # ==========================================
        # NEW: Backdoor Payload Variables
        # ==========================================
        self.is_backdoor = False
        self.backdoor_re = []
        self.backdoor_im = []

    def set_backdoor_rom(self, re_array, im_array):
        """Helper to cleanly pack a backdoor transaction."""
        self.is_backdoor = True
        self.backdoor_re = re_array
        self.backdoor_im = im_array

    def calculate_hw_registers(self):
        """Calculates generic equations based on the current f0, B, and cycles."""
        # Restored original math.floor logic
        self.FTW_start = int(math.floor((self.f0 * (1 << self.M)) / self.Fs))

        self.FTW_step = int(math.floor((self.B * (1 << self.M)) / (self.Fs * self.target_cycles)))

        # Mask to ensure they strictly fit within M bits
        bit_mask = (1 << self.M) - 1
        self.FTW_start &= bit_mask
        self.FTW_step  &= bit_mask

    def calculate_chirp(self, f0, B, cycles=4096):
        """Helper function to cleanly calculate a specific chirp without randomizing."""
        self.f0 = f0
        self.B = B
        self.target_cycles = cycles
        self.calculate_hw_registers()


    def set_bus_write(self, addr):
        """Maps the internal calculated values to the physical bus pins using memory addresses."""
        self.wr_en = 1
        self.rd_en = 0
        self.addr  = addr
        
        # Route the correct calculated data to the bus based on the target address
        if addr == 0x0:
            self.wr_data = self.FTW_start
        elif addr == 0x4:
            self.wr_data = self.FTW_step
        elif addr == 0x8:
            self.wr_data = self.target_cycles
        else:
            self.wr_data = 0  # Fallback for undefined addresses

    def set_bus_read(self, addr=0x0):
        """Sets up a read transaction to trigger the hardware dds_ready_flag."""
        self.wr_en = 0
        self.rd_en = 1
        self.addr  = addr
        self.wr_data = 0
        
    def set_bus_idle(self):
        """Clears the bus to prevent accidental reads or writes."""
        self.wr_en = 0
        self.rd_en = 0
        self.addr  = 0
        self.wr_data = 0

    def randomize(self):
        """Custom randomization for general testcases."""
        nyquist = self.Fs / 2
        self.f0 = random.uniform(0, nyquist / 2)             
        self.B  = random.uniform(1e6, nyquist - self.f0)     
        self.rst_n = random.choices([0, 1], weights=[5, 95])[0]

        # Trigger the math update using the newly randomized values
        self.calculate_hw_registers()

    def convert2string(self):
        return (f"{self.get_name()} "
                f"BUS[ addr={hex(self.addr)}, we={self.wr_en}, wdata={hex(self.wr_data)}, re={self.rd_en} ] -> "
                f"HW_STATE[ FTW_start={hex(self.FTW_start)}, FTW_step={hex(self.FTW_step)} ]")