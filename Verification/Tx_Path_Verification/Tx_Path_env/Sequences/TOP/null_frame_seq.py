import pyuvm
from pyuvm import *
from top_seq_item import top_item 
import cocotb

class null_frame_seq(uvm_sequence):
    """
    TC-001: Null Frame (Absolute Zero) Test.
    Validates that zero chirp (M_dds = 0) combined with zero OFDM data 
    produces a perfectly zeroed output, ensuring no DC offsets or noise.
    """
    def __init__(self, name="null_frame_seq"):
        super().__init__(name)
        self.num_frames = 1  
        
    async def body(self):
        
        # Target physical parameters for Zero Chirp
        f0_target = 0e6
        B_target  = 0e6  # Both 0 ensures FTW_start = 0 and FTW_step = 0
        
        # MAIN FRAME LOOP
        for i in range(self.num_frames):
            
            # ==========================================
            # 0. BACKDOOR MEMORY LOAD (Absolute Zero)
            # ==========================================
            # Guarantee the OFDM RAM is completely zeroed out
            req = top_item.create(f"req_backdoor_{i}")
            await self.start_item(req)
            
            # Create arrays of exactly 2048 zeros
            zeros_re = [0] * 2048 
            zeros_im = [0] * 2048
            
            req.set_backdoor_rom(zeros_re, zeros_im)
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
            # Cycle 6+: Bus goes idle, Hardware runs the pipeline with all zeros
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