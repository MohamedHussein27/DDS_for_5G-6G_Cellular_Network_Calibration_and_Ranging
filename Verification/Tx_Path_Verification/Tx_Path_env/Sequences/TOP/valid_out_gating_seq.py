import pyuvm
import random
from pyuvm import *
from top_seq_item import top_item
import cocotb

class valid_out_gating_seq(uvm_sequence):
    """
    TC-008: valid_out Gating Verification.
    Sends one frame. Observes valid_out timing relative to frame_start.
    Verifies valid_out is LOW between frames and during pipeline latency.
    Expected:
        - valid_out = 1 for exactly 4096 consecutive cycles per frame.
        - valid_out = 0 between frames.
    Pipeline latency L must remain constant (ASS-05 deterministic latency).
    Feature: F-04
    """
    def __init__(self, name="valid_out_gating_seq"):
        super().__init__(name)
        self.num_frames = 1  # Single frame is sufficient to verify gating

    async def body(self):

        # Nominal single-tone for clean valid_out observation (TC-008)
        f0_target = 50e6
        B_target  = 0e6   # Constant tone — isolate gating, not chirp behaviour

        # MAIN FRAME LOOP
        for i in range(self.num_frames):

            # ==========================================
            # 0. BACKDOOR MEMORY LOAD (Randomized)
            # ==========================================
            req = top_item.create(f"req_backdoor_{i}")
            await self.start_item(req)

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
            # 2. READ-BACK PHASE (Triggers Hardware)
            # ==========================================
            # Cycle 3: Read back START (addr 0x0)
            # valid_out must still be LOW here (pre-pipeline latency)
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
            # valid_out rises after pipeline latency L from this point
            req = top_item.create(f"req_r3_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x8)
            await self.finish_item(req)

            # ==========================================
            # 3. EXECUTION PHASE
            # Scoreboard must observe valid_out = 1 for exactly 4096 cycles
            # ==========================================
            for cycle in range(4096):
                req = top_item.create(f"req_exec_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)

            # ==========================================
            # 4. IDLE / INTER-FRAME GAP
            # Scoreboard must observe valid_out = 0 throughout this window
            # ==========================================
            for cycle in range(4096 * 5):
                req = top_item.create(f"req_idle_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)
