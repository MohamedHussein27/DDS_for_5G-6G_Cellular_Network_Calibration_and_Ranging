import pyuvm
from pyuvm import *
from top_seq_item import top_item 
import cocotb

class chirp_only_seq(uvm_sequence):
    """
    TC-002: Run 1 frame with a fixed chirp and zero OFDM symbols.
    """
    def __init__(self, name="chirp_only_seq"):
        super().__init__(name)
        self.num_frames = 1
        
    async def body(self):
        
        f0_target = 0e6
        B_target  = 200e6
        
        # MAIN FRAME LOOP
        for i in range(self.num_frames):
            
            # ==========================================
            # 0. BACKDOOR MEMORY LOAD (Via Driver)
            # ==========================================
            req = top_item.create(f"req_backdoor_{i}")
            await self.start_item(req)  # <-- Sequencer grants access first
            
            zeros_re = [0] * 2048 
            zeros_im = [0] * 2048
            req.set_backdoor_rom(zeros_re, zeros_im) # <-- Late assignment
            
            await self.finish_item(req) # <-- Send to driver

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
            # 2. READ-BACK PHASE (Triggers Hardware)
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

            # Cycle 5: Read back CYCLES (addr 0x8) -> Triggers dds_ready_flag
            req = top_item.create(f"req_r3_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x8)
            await self.finish_item(req)

            # ==========================================
            # 3. EXECUTION PHASE (Hardware Runs the Sweep)
            # ==========================================
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