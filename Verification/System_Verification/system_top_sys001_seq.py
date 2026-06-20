import pyuvm
from pyuvm import *
from system_top_seq_item import *
import cocotb
import random

class sys001_ofdm_loopback_Seq(uvm_sequence):
    """
    SYS-001: Single-Frame Loopback — OFDM.
    DDS chirp is held idle (f0 = 0 MHz, B = 0 MHz) so only the OFDM
    path is active. OFDM memory is backdoor-loaded with random 256-QAM
    I/Q symbols. RX is expected to recover the symbols through the
    equalizer with EVM within spec, and valid_out timing must align
    with the OFDM active window.
    """
    def __init__(self, name="sys001_ofdm_loopback_seq"):
        super().__init__(name)
        self.num_frames = 5

    async def body(self):

        # Target physical parameters — chirp disabled for this test
        f0_target = 0e6
        B_target  = 0e6

        # MAIN FRAME LOOP
        for i in range(self.num_frames):

            # ==========================================
            # 0. BACKDOOR MEMORY LOAD (Randomized OFDM symbols)
            # ==========================================
            req = system_top_item.create(f"req_backdoor_{i}")
            await self.start_item(req)

            # Generate 2048 random 16-bit signed integers (256-QAM I/Q)
            rand_re = [random.randint(-32768, 32767) for _ in range(2048)]
            rand_im = [random.randint(-32768, 32767) for _ in range(2048)]

            req.set_backdoor_rom(rand_re, rand_im)
            await self.finish_item(req)

            # ==========================================
            # 1. WRITE PHASE (Configuration)
            # ==========================================
            # Cycle 0: Write FTW_START (addr 0x0)
            req = system_top_item.create(f"req_w1_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_write(addr=0x0)
            await self.finish_item(req)

            # Cycle 1: Write FTW_STEP (addr 0x4)
            req = system_top_item.create(f"req_w2_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_write(addr=0x4)
            await self.finish_item(req)

            # Cycle 2: Write CYCLES (addr 0x8)
            req = system_top_item.create(f"req_w3_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_write(addr=0x8)
            await self.finish_item(req)

            # ==========================================
            # 2. READ-BACK PHASE (Verification & Trigger)
            # ==========================================
            # Cycle 3: Read back START (addr 0x0)
            req = system_top_item.create(f"req_r1_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x0)
            await self.finish_item(req)

            # Cycle 4: Read back STEP (addr 0x4)
            req = system_top_item.create(f"req_r2_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x4)
            await self.finish_item(req)

            # Cycle 5: Read back CYCLES (addr 0x8) -> Triggers hardware
            req = system_top_item.create(f"req_r3_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x8)
            await self.finish_item(req)

            # ==========================================
            # 3. EXECUTION PHASE (Active Data Cycles)
            # ==========================================
            # Hardware runs the pipeline: TX OFDM frame -> channel -> RX recovery
            for cycle in range(4096):
                req = system_top_item.create(f"req_exec_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)

            # ==========================================
            # 4. IDLE / INTER-FRAME GAP
            # ==========================================
            for cycle in range(4096 * 5):
                req = system_top_item.create(f"req_idle_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)
