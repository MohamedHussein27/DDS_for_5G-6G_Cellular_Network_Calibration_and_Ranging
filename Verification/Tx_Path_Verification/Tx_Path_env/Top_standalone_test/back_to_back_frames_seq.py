import pyuvm
import random
from pyuvm import *
from top_seq_item import top_item
import cocotb

class back_to_back_frames_seq(uvm_sequence):
    """
    TC-012: Back-to-Back 20 Frames — No Gaps. [HIGH Priority]
    Transmit 20 consecutive frames with zero idle cycles between them.
    Each frame uses a different random OFDM seed.
    Expected: For every one of the 20 frames, no inter-frame contamination
    observable at tx_out.
    PASS: 0 errors across all 20 frames compared to the fixed-point Python model.
    """
    def __init__(self, name="back_to_back_frames_seq"):
        super().__init__(name)
        self.num_frames = 20

    async def body(self):

        # Nyquist boundary for safe random f0 selection
        Fs      = 491520000
        nyquist = Fs / 2

        # MAIN FRAME LOOP — 20 frames, zero idle gap between them
        for i in range(self.num_frames):

            # Each frame gets a unique random OFDM seed / chirp parameters
            f0_target = random.uniform(10e6, nyquist / 2)
            B_target  = random.uniform(1e6,  nyquist / 2 - f0_target + 1e6)

            # ==========================================
            # 0. BACKDOOR MEMORY LOAD (Unique random seed per frame)
            # ==========================================
            req = top_item.create(f"req_backdoor_{i}")
            await self.start_item(req)

            # Different random OFDM symbols every frame — stresses inter-frame isolation
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
            # 3. EXECUTION PHASE
            # Zero idle cycles between frames — next frame config starts immediately
            # ==========================================
            for cycle in range(4096):
                req = top_item.create(f"req_exec_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)

        # ==========================================
        # 4. IDLE / INTER-FRAME GAP just to get the output for the last frame
        # ==========================================
        for cycle in range(4096 * 5):
            req = top_item.create(f"req_idle_{i}_{cycle}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_idle()
            await self.finish_item(req)  