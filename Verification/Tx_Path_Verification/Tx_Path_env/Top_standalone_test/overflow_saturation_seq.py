import pyuvm
import random
from pyuvm import *
from top_seq_item import top_item
import cocotb

class overflow_saturation_seq(uvm_sequence):
    """
    TC-007: Overflow / Saturation Behaviour.
    Drives worst-case maximum-amplitude stimulus: max M_dds + all OFDM corner
    symbols simultaneously.
    Expected: If any sample hits ±32767, the next sample must NOT jump to the
    opposite extreme (no wrap-around). Saturate-don't-wrap.
    Feature: F-02
    """
    def __init__(self, name="overflow_saturation_seq"):
        super().__init__(name)
        self.num_frames = 1

    async def body(self):

        # Worst-case maximum-amplitude: frequency near Nyquist to stress amplitude path
        f0_target = 122880000e0  # ~Fs/4, maximises DDS output amplitude
        B_target  = 0e6          # Constant tone — isolate saturation, not chirp

        # MAIN FRAME LOOP
        for i in range(self.num_frames):

            # ==========================================
            # 0. BACKDOOR MEMORY LOAD (Full-Scale Corner Symbols)
            # ==========================================
            req = top_item.create(f"req_backdoor_{i}")
            await self.start_item(req)

            # Worst-case OFDM corner: all symbols at full scale ±32767
            max_re = [32767  if j % 2 == 0 else -32767 for j in range(2048)]
            max_im = [-32767 if j % 2 == 0 else  32767 for j in range(2048)]

            req.set_backdoor_rom(max_re, max_im)
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
            # 3. EXECUTION PHASE (Hardware Runs — Monitor for Wrap-Around)
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
