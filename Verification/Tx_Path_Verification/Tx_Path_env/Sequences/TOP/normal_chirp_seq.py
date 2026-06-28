"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: normal_chirp_seq.py

    Description:
        Sequence that runs a single linear chirp frame with:
            f0 = 0 Hz
            B  = 200 MHz
        No backdoor memory access – the OFDM memory is used exactly as it is.
        Validates the full DDS TX datapath without random data injection.
"""
import pyuvm
from pyuvm import *
from top_seq_item import top_item

class normal_chirp_seq(uvm_sequence):
    """
    TC-005: Run 1 frame with the default 200 MHz linear chirp (f0=0).
    Uses the OFDM memory content already present – no backdoor loading.
    """
    def __init__(self, name="normal_chirp_seq"):
        super().__init__(name)
        self.num_frames = 1          # Change to 100 if desired

    async def body(self):
        # Physical parameters for the chirp
        f0_target = 0.0
        B_target  = 200e6

        for i in range(self.num_frames):
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