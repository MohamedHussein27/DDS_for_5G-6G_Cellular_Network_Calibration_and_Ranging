import pyuvm
import random
from pyuvm import *
from top_seq_item import top_item
import cocotb

class long_gap_between_frames_seq(uvm_sequence):
    """
    TC-013: Long Gap Between Frames. [MED Priority]
    Send frame 1. Insert 20,000 idle cycles (no frame_start). Send frame 2
    with a different stimulus.
    Expected: Frame 2 has no stale pipeline state from frame 1 corrupting it.
    PASS: Frame 2 outputs exactly match the Python reference (Error = 0).
    """
    def __init__(self, name="long_gap_between_frames_seq"):
        super().__init__(name)
        self.long_gap_cycles = 20000  # TC-013 specification

    async def body(self):

        Fs      = 491520000
        nyquist = Fs / 2

        # ============================================================
        # FRAME 1 — Arbitrary stimulus
        # ============================================================
        f0_frame1 = 40e6
        B_frame1  = 20e6

        # 0. Backdoor load for frame 1
        req = top_item.create("req_backdoor_f1")
        await self.start_item(req)
        # Generate 2048 random 16-bit signed integers for Real and Imaginary.
        # (Matches your active subcarrier depth. Covers the 1500 symbols 
        # mentioned in TC-003 while safely initializing the entire array).
        # max real and max imag in conestellation: 1.1504 -> q8.8 -> 295
        rand_re_1 = [random.randint(-295, 295) for _ in range(2048)]
        rand_im_1 = [random.randint(-295, 295) for _ in range(2048)]
        #req.set_backdoor_rom(rand_re_1, rand_im_1)
        await self.finish_item(req)

        # 1. Write phase — frame 1
        req = top_item.create("req_f1_w1")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame1, B=B_frame1)
        req.set_bus_write(addr=0x0)
        await self.finish_item(req)

        req = top_item.create("req_f1_w2")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame1, B=B_frame1)
        req.set_bus_write(addr=0x4)
        await self.finish_item(req)

        req = top_item.create("req_f1_w3")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame1, B=B_frame1)
        req.set_bus_write(addr=0x8)
        await self.finish_item(req)

        # 2. Read-back phase — frame 1 (triggers dds_ready_flag)
        req = top_item.create("req_f1_r1")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame1, B=B_frame1)
        req.set_bus_read(addr=0x0)
        await self.finish_item(req)

        req = top_item.create("req_f1_r2")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame1, B=B_frame1)
        req.set_bus_read(addr=0x4)
        await self.finish_item(req)

        req = top_item.create("req_f1_r3")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame1, B=B_frame1)
        req.set_bus_read(addr=0x8)
        await self.finish_item(req)

        # 3. Execution phase — frame 1
        for cycle in range(4096):
            req = top_item.create(f"req_f1_exec_{cycle}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_frame1, B=B_frame1)
            req.set_bus_idle()
            await self.finish_item(req)

        # ============================================================
        # LONG GAP — 20,000 idle cycles, no frame_start asserted
        # Pipeline must drain and hold no stale state during this window
        # ============================================================
        for cycle in range(self.long_gap_cycles):
            req = top_item.create(f"req_long_gap_{cycle}")
            await self.start_item(req)
            req.set_bus_idle()
            await self.finish_item(req)

        # ============================================================
        # FRAME 2 — Different stimulus (stale-state stress target)
        # ============================================================
        f0_frame2 = 90e6   # Deliberately different from frame 1
        B_frame2  = 15e6

        # 0. Backdoor load for frame 2 — completely different OFDM seed
        req = top_item.create("req_backdoor_f2")
        await self.start_item(req)
        rand_re_f2 = [random.randint(-295, 295) for _ in range(2048)]
        rand_im_f2 = [random.randint(-295, 295) for _ in range(2048)]
        #req.set_backdoor_rom(rand_re_f2, rand_im_f2)
        await self.finish_item(req)

        # 1. Write phase — frame 2
        req = top_item.create("req_f2_w1")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame2, B=B_frame2)
        req.set_bus_write(addr=0x0)
        await self.finish_item(req)

        req = top_item.create("req_f2_w2")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame2, B=B_frame2)
        req.set_bus_write(addr=0x4)
        await self.finish_item(req)

        req = top_item.create("req_f2_w3")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame2, B=B_frame2)
        req.set_bus_write(addr=0x8)
        await self.finish_item(req)

        # 2. Read-back phase — frame 2 (triggers dds_ready_flag)
        req = top_item.create("req_f2_r1")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame2, B=B_frame2)
        req.set_bus_read(addr=0x0)
        await self.finish_item(req)

        req = top_item.create("req_f2_r2")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame2, B=B_frame2)
        req.set_bus_read(addr=0x4)
        await self.finish_item(req)

        req = top_item.create("req_f2_r3")
        await self.start_item(req)
        req.calculate_chirp(f0=f0_frame2, B=B_frame2)
        req.set_bus_read(addr=0x8)
        await self.finish_item(req)

        # 3. Execution phase — frame 2 (scoreboard compares against Python reference)
        for cycle in range(4096):
            req = top_item.create(f"req_f2_exec_{cycle}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_frame2, B=B_frame2)
            req.set_bus_idle()
            await self.finish_item(req)

        # 4. Final idle gap
        for cycle in range(4096 * 5):
            req = top_item.create(f"req_final_idle_{cycle}")
            await self.start_item(req)
            req.set_bus_idle()
            await self.finish_item(req)
