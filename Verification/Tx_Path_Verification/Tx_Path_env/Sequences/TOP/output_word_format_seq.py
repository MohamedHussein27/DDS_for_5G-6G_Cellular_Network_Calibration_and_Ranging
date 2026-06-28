import pyuvm
import random
from pyuvm import *
from top_seq_item import top_item
import cocotb

class output_word_format_seq(uvm_sequence):
    """
    TC-006: Output Word Format — Normal Operation.
    Dual-band nominal stimulus. Inspects all 4096 tx_out samples for correct
    signed 16-bit fixed-point format (WL=16, FL=5).
    All samples must be within [-32768, +32767].
    No sample equals -32768 unless saturation is explicitly expected.
    Feature: F-02
    """
    def __init__(self, name="output_word_format_seq"):
        super().__init__(name)
        self.num_frames = 1

    async def body(self):

        # Dual-band nominal stimulus parameters (TC-006)
        # Band 1: Low-frequency tone, Band 2: Mid-frequency tone
        f0_band1 = 30e6
        f0_band2 = 80e6
        B_target  = 0e6   # Constant tones for format inspection

        # MAIN FRAME LOOP
        for i in range(self.num_frames):

            # ==========================================
            # 0. BACKDOOR MEMORY LOAD (Nominal Values)
            # ==========================================
            req = top_item.create(f"req_backdoor_{i}")
            await self.start_item(req)

            # Nominal 16-bit signed values — not full-scale, not zero
            nom_re = [random.randint(-16384, 16383) for _ in range(2048)]
            nom_im = [random.randint(-16384, 16383) for _ in range(2048)]

            req.set_backdoor_rom(nom_re, nom_im)
            await self.finish_item(req)

            # ==========================================
            # BAND 1 FRAME
            # ==========================================
            # Cycle 0: Write FTW_START (addr 0x0)
            req = top_item.create(f"req_b1_w1_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band1, B=B_target)
            req.set_bus_write(addr=0x0)
            await self.finish_item(req)

            # Cycle 1: Write FTW_STEP (addr 0x4)
            req = top_item.create(f"req_b1_w2_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band1, B=B_target)
            req.set_bus_write(addr=0x4)
            await self.finish_item(req)

            # Cycle 2: Write CYCLES (addr 0x8)
            req = top_item.create(f"req_b1_w3_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band1, B=B_target)
            req.set_bus_write(addr=0x8)
            await self.finish_item(req)

            # Cycle 3: Read back START (addr 0x0)
            req = top_item.create(f"req_b1_r1_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band1, B=B_target)
            req.set_bus_read(addr=0x0)
            await self.finish_item(req)

            # Cycle 4: Read back STEP (addr 0x4)
            req = top_item.create(f"req_b1_r2_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band1, B=B_target)
            req.set_bus_read(addr=0x4)
            await self.finish_item(req)

            # Cycle 5: Read back CYCLES (addr 0x8) -> Triggers dds_ready_flag
            req = top_item.create(f"req_b1_r3_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band1, B=B_target)
            req.set_bus_read(addr=0x8)
            await self.finish_item(req)

            # Execution Phase — observe all 4096 output samples
            for cycle in range(4096):
                req = top_item.create(f"req_b1_exec_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_band1, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)

            # Inter-frame gap
            for cycle in range(4096 * 2):
                req = top_item.create(f"req_b1_idle_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_band1, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)

            # ==========================================
            # BAND 2 FRAME
            # ==========================================
            # Cycle 0: Write FTW_START (addr 0x0)
            req = top_item.create(f"req_b2_w1_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band2, B=B_target)
            req.set_bus_write(addr=0x0)
            await self.finish_item(req)

            # Cycle 1: Write FTW_STEP (addr 0x4)
            req = top_item.create(f"req_b2_w2_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band2, B=B_target)
            req.set_bus_write(addr=0x4)
            await self.finish_item(req)

            # Cycle 2: Write CYCLES (addr 0x8)
            req = top_item.create(f"req_b2_w3_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band2, B=B_target)
            req.set_bus_write(addr=0x8)
            await self.finish_item(req)

            # Cycle 3: Read back START (addr 0x0)
            req = top_item.create(f"req_b2_r1_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band2, B=B_target)
            req.set_bus_read(addr=0x0)
            await self.finish_item(req)

            # Cycle 4: Read back STEP (addr 0x4)
            req = top_item.create(f"req_b2_r2_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band2, B=B_target)
            req.set_bus_read(addr=0x4)
            await self.finish_item(req)

            # Cycle 5: Read back CYCLES (addr 0x8) -> Triggers dds_ready_flag
            req = top_item.create(f"req_b2_r3_{i}")
            await self.start_item(req)
            req.calculate_chirp(f0=f0_band2, B=B_target)
            req.set_bus_read(addr=0x8)
            await self.finish_item(req)

            # Execution Phase — observe all 4096 output samples
            for cycle in range(4096):
                req = top_item.create(f"req_b2_exec_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_band2, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)

            # Inter-frame gap
            for cycle in range(4096 * 5):
                req = top_item.create(f"req_b2_idle_{i}_{cycle}")
                await self.start_item(req)
                req.calculate_chirp(f0=f0_band2, B=B_target)
                req.set_bus_idle()
                await self.finish_item(req)
