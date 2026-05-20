import pyuvm
from pyuvm import *
from top_seq_item import top_item 
import cocotb
import random

class full_sweep_seq(uvm_sequence):
    """
    TC-005: Full M_dds Sweep.
    Applies a sweeping chirp (f0 = 10 MHz, B = 200 MHz).
    The OFDM memory is fully randomized to stress the mixed-signal datapath.
    """
    def __init__(self, name="full_sweep_seq"):
        super().__init__(name)
        self.num_frames = 1  
        
    async def body(self):
        
        # Target physical parameters
        f0_target = 10e6
        B_target  = 200e6
        
        # MAIN FRAME LOOP
        for i in range(self.num_frames):
            
            # ==========================================
            # 0. BACKDOOR MEMORY LOAD (Randomized)
            # ==========================================
            req = top_item.create(f"req_backdoor_{i}")
            await self.start_item(req)
            
            # Generate 2048 random 16-bit signed integers
            rand_re = [random.randint(-32768, 32767) for _ in range(2048)]
            rand_im = [random.randint(-32768, 32767) for _ in range(2048)]
            
            req.set_backdoor_rom(rand_re, rand_im)
            await self.finish_item(req)

            # ==========================================
            # 1. WRITE PHASE (Configuration)
            # ==========================================
            # Cycle 0: Write FTW_START (addr 0x0)
            req = top_item.create(f"req_w1_{i}")
            await self.start_item(req) 
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_write(addr=0x0)
            await self.finish_item(req)

            # Cycle 1: Write FTW_STEP (addr 0x4)
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

            # Cycle 5: Read back CYCLES (addr 0x8) -> Triggers hardware
            req = top_item.create(f"req_r3_{i}")
            await self.start_item(req) 
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x8)
            await self.finish_item(req)

            # ==========================================
            # 3. EXECUTION PHASE (Active Data Cycles)
            # ==========================================
            # Hardware runs the pipeline with the sweep and random OFDM
            for cycle in range(4096): 
                req = top_item.create(f"req_exec_{i}_{cycle}")
                await self.start_item(req) 
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)
            
            # ==========================================
            # 4. IDLE / INTER-FRAME GAP
            # ==========================================
            for cycle in range(4096 * 5):
                req = top_item.create(f"req_idle_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)