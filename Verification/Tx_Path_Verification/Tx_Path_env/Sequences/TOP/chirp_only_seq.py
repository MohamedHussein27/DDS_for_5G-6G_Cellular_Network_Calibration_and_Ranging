import pyuvm
from pyuvm import *
from top_seq_item import top_item 
import cocotb

class chirp_only_seq(uvm_sequence):
    """
    TC-002: Run 1 frame with a sweeping chirp and zero OFDM symbols.
    This acts as the ultimate stress test for the TX Top Datapath.
    """
    def __init__(self, name="chirp_only_seq"):
        super().__init__(name)
        self.num_frames = 1
        
    def backdoor_zero_rom(self):
        """
        Bypasses the physical bus to dynamically clear the internal OFDM ROM.
        Forces all OFDM symbols to 0 for a pure chirp test (TC-002).
        """
        # 1. Fetch the DUT from the UVM config database
        dut = ConfigDB().get(self, "", "DUT")
        

        dut_ram_re = self.dut.u_ofdm_rom.rom_real 
        dut_ram_im = self.dut.u_ofdm_rom.rom_imag 
        
        ram_depth = len(dut_ram_re)
        
        # 2. Loop through and force all values to strictly 0
        for i in range(ram_depth):
            dut_ram_re[i].value = 0
            dut_ram_im[i].value = 0
            
        uvm_root().logger.info(f"BACKDOOR: Successfully zeroed out {ram_depth} entries in OFDM ROM.")

    async def body(self):
        
        # Target physical parameters
        f0_target = 30e6
        B_target  = 200e6
        
        # MAIN FRAME LOOP
        for i in range(self.num_frames):
            
            # ==========================================
            # 0. BACKDOOR MEMORY LOAD
            # ==========================================
            # Guarantee a completely blank OFDM map before we trigger the sweep
            self.backdoor_zero_rom()
            
            # ==========================================
            # 1. WRITE PHASE (Configuration)
            # ==========================================
            # Cycle 0: Write FTW_START (addr 0x0)
            req = top_item.create(f"req_w1_{i}")
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_write(addr=0x0)
            await self.start_item(req); await self.finish_item(req)

            # Cycle 1: Write FTW_STEP (addr 0x4)
            req = top_item.create(f"req_w2_{i}")
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_write(addr=0x4)
            await self.start_item(req); await self.finish_item(req)

            # Cycle 2: Write CYCLES (addr 0x8)
            req = top_item.create(f"req_w3_{i}")
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_write(addr=0x8)
            await self.start_item(req); await self.finish_item(req)

            # ==========================================
            # 2. READ-BACK PHASE (Triggers Hardware)
            # ==========================================
            # Cycle 3: Read back START (addr 0x0)
            req = top_item.create(f"req_r1_{i}")
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x0)
            await self.start_item(req); await self.finish_item(req)

            # Cycle 4: Read back STEP (addr 0x4)
            req = top_item.create(f"req_r2_{i}")
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x4)
            await self.start_item(req); await self.finish_item(req)

            # Cycle 5: Read back CYCLES (addr 0x8) -> Triggers dds_ready_flag
            req = top_item.create(f"req_r3_{i}")
            req.calculate_chirp(f0=f0_target, B=B_target)
            req.set_bus_read(addr=0x8)
            await self.start_item(req); await self.finish_item(req)

            # ==========================================
            # 3. EXECUTION PHASE (Hardware Runs the Sweep)
            # ==========================================
            # Cycle 6+: Bus goes idle, Hardware runs the sweep!
            for cycle in range(4096):
                req = top_item.create(f"req_exec_{i}_{cycle}")
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.start_item(req); await self.finish_item(req)
            
            # ==========================================
            # 4. IDLE / INTER-FRAME GAP
            # ==========================================
            # Stable idle state between frames
            for cycle in range(4096 * 5):
                req = top_item.create(f"req_idle_{i}_{cycle}")
                req.calculate_chirp(f0=f0_target, B=B_target)
                req.set_bus_idle()
                await self.start_item(req); await self.finish_item(req)