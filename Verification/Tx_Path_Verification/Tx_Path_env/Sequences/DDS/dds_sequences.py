import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 

class dds_reset(uvm_sequence):
    def __init__(self, name="dds_reset"):
        super().__init__(name)

    async def body(self):
        req = dds_seq_item("req")
        await self.start_item(req)
        req.randomize()
        req.rst_n = 0  # FORCE RESET ACTIVE
        await self.finish_item(req)                                                  
                                    

class dds_random_seq(dds_reset):
    def __init__(self, name="dds_random_seq"):
        super().__init__(name)
        # Default to 10, but allow the Test to overwrite this!
        self.num_items = 10 
        
    async def body(self):
        await super.body()  # Ensure we start with a reset transaction to initialize the hardware
        for _ in range(self.num_items):
            req = dds_seq_item("req")
            await self.start_item(req)
            req.randomize()
            await self.finish_item(req)
 
 
class dds_singletone_seq(dds_reset):
    def __init__(self, name="dds_singletone_seq"):
        super().__init__(name)
        self.num_items = 15
    async def body(self):
        await super.body()  # Start with a reset to ensure clean state
        # Generate 15 pure single tones (No chirp acceleration)
        for _ in range(self.num_items):
            req = dds_seq_item("req")
            await self.start_item(req)
            # Force FTW_step to 0, randomize the rest
            req.randomize()
            req.rst_n = 1
            req.FTW_step = 0
            await self.finish_item(req)
            
class dds_cyclestress_seq(dds_reset):
    def __init__(self, name="dds_cyclestress_seq"):
        super().__init__(name)
    async def body(self):
        await super.body()
        # Force the cycle counter to overflow almost instantly
        short_cycles = [1, 2, 3, 5]
        for c in short_cycles:
            req = dds_seq_item("req")
            await self.start_item(req)
            req.randomize()
            req.rst_n = 1
            req.cycles = c  # Force short cycle counts to stress the counter logic
            await self.finish_item(req)                    

class dds_chirpsweep_seq(dds_reset):
    
    def __init__(self, name="dds_chirpsweep_seq"):
        super().__init__(name)
    async def body(self):
        await super.body()
        # Array of aggressive tuning steps
        chirp_steps = [0x100, 0x5000, 0x1FFFF, 0x8FFFF]
        for step_val in chirp_steps:
            req = dds_seq_item("req")
            await self.start_item(req)
            # Fix the cycles so the wave has time to propagate, force the step
            req.randomize()
            req.rst_n = 1
            req.cycles = 1000  # Give enough time for the chirp to evolve
            req.FTW_step = step_val  # Force aggressive chirp steps
            await self.finish_item(req) 
            
class dds_reset_recovery_seq(dds_reset):
    def __init__(self, name="dds_reset_recovery_seq"):
        super().__init__(name)

    async def body(self):
        await super.body()
        for _ in range(5):
            # Step 1: Send a reset state to kill whatever the hardware is doing
            req_kill = dds_seq_item("req_kill")
            await self.start_item(req_kill)
            req_kill.randomize()
            req_kill.rst_n = 0  # FORCE RESET ACTIVE
            await self.finish_item(req_kill)
            
            # Step 2: Immediately send a valid transaction to ensure it woke up
            req_recover = dds_seq_item("req_recover")
            await self.start_item(req_recover)
            req_recover.randomize()
            req_recover.rst_n = 1 # FORCE RESET INACTIVE
            await self.finish_item(req_recover)
                                   
class dds_cornercase_seq(dds_reset):
    def __init__(self, name="dds_cornercase_seq"):
        super().__init__(name)

    async def body(self):
        await super.body()
        # Test Absolute Maximums 
        req = dds_seq_item("req")
        await self.start_item(req)
        req.randomize()
        
        req.rst_n = 1
        req.FTW_start = 0xFFFFFFFF  # All 1s
        req.FTW_step = 0xFFFFFFFF   # All 1s
        req.cycles = 0x1FFF         # Max 13-bit value (8191)
        
        await self.finish_item(req)
        
        
class dds_fft_boundary_seq(dds_reset):
    def __init__(self, name="dds_fft_boundary_seq"):
        super().__init__(name)

    async def body(self):
        await super.body()
        # 4095 = 1 cycle short (Underflow)
        # 4096 = Perfect FFT Symbol
        # 4097 = 1 cycle over (Overflow/Backpressure)
        boundaries = [4095, 4096, 4097]
        
        for c in boundaries:
            req = dds_seq_item("req")
            await self.start_item(req)
            req.randomize()
            
            req.rst_n = 1
            req.cycles = c  # Force the specific FFT boundary
            
            await self.finish_item(req)        
            
            
   