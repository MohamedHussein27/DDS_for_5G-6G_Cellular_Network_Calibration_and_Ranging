import pyuvm
from pyuvm import *
from top_seq_item import top_item 
import cocotb
import random

class ofdm_only_seq(uvm_sequence):
    """
    TC-003: OFDM-Only — Zero Chirp.
    Validates the TX datapath when the DDS is completely silent (M_dds = 0),
    ensuring only the randomized OFDM symbols propagate through the pipeline.
    """
    def __init__(self, name="ofdm_only_seq"):
        super().__init__(name)
        self.num_frames = 1  
        
    async def body(self):
        
        # Target physical parameters for Zero Chirp (TC-003)
        f0_target = 0e6
        B_target  = 0e6  # Both 0 ensures FTW_start = 0 and FTW_step = 0
        
        # MAIN FRAME LOOP
        for i in range(self.num_frames):
            
            # ==========================================
            # 0. BACKDOOR MEMORY LOAD (Randomized!)
            # ==========================================
            # Guarantee the OFDM RAM is filled with random data before triggering
            req = top_item.create(f"req_backdoor_{i}")
            await self.start_item(req)
            
            # Generate 2048 random 16-bit signed integers for Real and Imaginary.
            # (Matches your active subcarrier depth. Covers the 1500 symbols 
            # mentioned in TC-003 while safely initializing the entire array).
            # max real and max imag in conestellation: 1.1504 -> q8.8 -> 295
            rand_re = [random.randint(-295, 295) for _ in range(2048)]
            rand_im = [random.randint(-295, 295) for _ in range(2048)]
            
            req.set_backdoor_rom(rand_re, rand_im)
            await self.finish_item(req)

            # ==========================================
            # 1. WRITE PHASE (Configuration)
            # ==========================================
            # Cycle 0: Write FTW_START (addr 0x0) -> Will be 0
            req = top_item.create(f"req_w1_{i}")
            await self.start_item(req) 
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_write(addr=0x0)
            await self.finish_item(req)

            # Cycle 1: Write FTW_STEP (addr 0x4) -> Will be 0
            req = top_item.create(f"req_w2_{i}")
            await self.start_item(req) 
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_write(addr=0x4)
            await self.finish_item(req)

            # Cycle 2: Write CYCLES (addr 0x8)
            req = top_item.create(f"req_w3_{i}")
            await self.start_item(req) 
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_write(addr=0x8)
            await self.finish_item(req)

            # ==========================================
            # 2. READ-BACK PHASE (Verification & Trigger)
            # ==========================================
            # Cycle 3: Read back START (addr 0x0)
            req = top_item.create(f"req_r1_{i}")
            await self.start_item(req) 
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x0)
            await self.finish_item(req)

            # Cycle 4: Read back STEP (addr 0x4)
            req = top_item.create(f"req_r2_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x4)
            await self.finish_item(req)

            # Cycle 5: Read back CYCLES (addr 0x8)
            # This 3rd read triggers 'dds_ready_flag' inside the hardware!
            req = top_item.create(f"req_r3_{i}")
            await self.start_item(req) 
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x8)
            await self.finish_item(req)

            # ==========================================
            # 3. EXECUTION PHASE (Active Data Cycles)
            # ==========================================
            # Cycle 6+: Bus goes idle, Hardware runs the pipeline with M_dds = 0
            for cycle in range(4096): 
                req = top_item.create(f"req_exec_{i}_{cycle}")
                await self.start_item(req) 
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)
            
            # ==========================================
            # 4. IDLE / INTER-FRAME GAP
            # ==========================================
            # Hold the bus idle between frames
            for cycle in range(4096 * 5):
                req = top_item.create(f"req_idle_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)