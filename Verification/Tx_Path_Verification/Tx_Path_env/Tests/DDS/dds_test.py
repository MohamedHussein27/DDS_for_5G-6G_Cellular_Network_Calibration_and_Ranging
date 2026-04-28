"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

Module: DDS Full Regression Test Suite
Type:   pyuvm Test Class
Target: DDS Phase Accumulator (5G ISAC Transmitter)

Description:
    This is the top-level UVM regression test for the Direct Digital Synthesizer 
    (DDS) block. It is designed to hit 100% functional coverage by exercising 
    the hardware through both directed edge-cases and high-volume random traffic.

Methodology ("Fail-Fast" Execution):
    The run_phase executes sequences in an escalating order of complexity. 
    It intentionally targets known hardware weak points (resets, overflows, 
    and boundaries) before attempting continuous mathematical sweeps. This 
    guarantees that fundamental digital logic bugs are caught immediately 
    before wasting simulation time on high-volume random regression.

Test Sequence Flow:
    1. Reset Recovery   : Intentionally yanks reset low during operation and 
                          verifies the accumulator clears gracefully.
    2. Corner Cases     : Forces 32-bit registers to maximum limits (0xFFFFFFFF) 
                          to check adder/carry-out boundaries.
    3. FFT Boundaries   : Tests execution lengths immediately around the 4096-point 
                          FFT boundary (4095, 4096, 4097) to stress FIFOs.
    4. Single Tone (LO) : Verifies pure LO frequency generation (FTW_step = 0).
    5. Cycle Stress     : Spams the hardware with micro-transactions (cycles < 5) 
                          to stress the internal cycle counter state machine.
    6. Chirp Sweep      : Simulates 5G ISAC radar sweeps using highly aggressive 
                          frequency steps to verify quadrant folding logic.
    7. Random Volume    : High-volume, constrained-random stimulus to hunt for 
                          hidden topological bugs.
    8. Full Regression     : Executes all of the above in a single run for final verification sign-off.                      


"""
import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

# Import the base test and the specific sequences
from base_test import base_test
from dds_sequences import * 


# ---------------------------------------------------------
# Test 1: The Full Regression (Runs everything)
# ---------------------------------------------------------
@pyuvm.test()
class dds_full_regression_test(base_test):
    def build_phase(self):
        super().build_phase()
        ConfigDB.set(self, "*", "VERIF_MODE", "DDS")  # Set the verification mode to TOP for this test
        self.reset_seq = dds_reset_recovery_seq("reset_seq")
        self.random_seq = dds_random_seq("random_seq")
        
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        await self.generate_clock()
        await self.run_initial_setup()
        self.logger.info(f"Starting sequence: {self.reset_seq.get_name()}")
        await self.reset_seq.start(self.env.agt.sqr)
        self.logger.info(f"Finished sequence: {self.reset_seq.get_name()}")
        self.logger.info(f"Starting sequence: {self.random_seq.get_name()}")
        await self.random_seq.start(self.env.agt.sqr)
        self.logger.info(f"Finished sequence: {self.random_seq.get_name()}")
        self.drop_objection()

# ---------------------------------------------------------
# Test 2: Just the ISAC Chirp Sweep (For fast debugging)
# ---------------------------------------------------------
@pyuvm.test()
class dds_chirp_only_test(base_test):
    def build_phase(self):
        super().build_phase()
        self.chirp_seq = dds_chirpsweep_seq("chirp_seq")
        
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        await self.generate_clock()
        await self.run_initial_setup()
        self.logger.info(f"Starting sequence: {self.chirp_seq.get_name()}")
        # Runs ONLY the chirp sequence
        await self.chirp_seq.start(self.env.agt.sqr)
        self.logger.info(f"finished sequence: {self.chirp_seq.get_name()}")
        self.drop_objection() 
        
# ---------------------------------------------------------
# Test 3: cyclestress test (For fast debugging of cycle counter overflow handling)
# ---------------------------------------------------------
@pyuvm.test()
class dds_cyclestress_test(base_test):
    def build_phase(self):
        super().build_phase()
        self.cyclestress_seq = dds_cyclestress_seq("cyclestress_seq")
        
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        await self.generate_clock()
        await self.run_initial_setup()
        self.logger.info(f"Starting sequence: {self.cyclestress_seq.get_name()}")
        # Runs ONLY the chirp sequence
        await self.cyclestress_seq.start(self.env.agt.sqr)
        self.logger.info(f"finished sequence: {self.cyclestress_seq.get_name()}")
        self.drop_objection() 
# ---------------------------------------------------------
# Test 4: Singletone only test (For fast debugging of single tone generation and FTW_step handling)
# ---------------------------------------------------------        
@pyuvm.test()
class dds_singletone_only_test(base_test):
    def build_phase(self):
        super().build_phase()
        self.singletone_seq = dds_singletone_seq("singletone_seq")
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        await self.generate_clock()
        await self.run_initial_setup()
        self.logger.info(f"Starting sequence: {self.singletone_seq.get_name()}")
        # Runs ONLY the chirp sequence
        await self.singletone_seq.start(self.env.agt.sqr)
        self.logger.info(f"finished sequence: {self.singletone_seq.get_name()}")
        self.drop_objection() 
 # ---------------------------------------------------------
# Test 5: Reset Recovery only test (For fast debugging of reset behavior and recovery)
# ---------------------------------------------------------        

@pyuvm.test()
class dds_resetrecovery_only_test(base_test):
    def build_phase(self):
        super().build_phase()
        self.reset_recovery_seq = dds_reset_recovery_seq.create("reset_recovery_seq")
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        await self.generate_clock()
        await self.run_initial_setup()
        self.logger.info(f"Starting sequence: {self.reset_recovery_seq.get_name()}")
        # Runs ONLY the chirp sequence
        await self.reset_recovery_seq.start(self.env.agt.sqr)
        self.logger.info(f"finished sequence: {self.reset_recovery_seq.get_name()}")
        self.drop_objection()   
# ---------------------------------------------------------
# Test 6: Corner Case only test (For fast debugging of corner cases and edge conditions)
# ---------------------------------------------------------        
@pyuvm.test()
class dds_cornercase_only_test(base_test):
    def build_phase(self):
        super().build_phase()
        self.cornercase_seq= dds_cornercase_seq.create("cornercase_seq")
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        await self.generate_clock()
        await self.run_initial_setup()
        self.logger.info(f"Starting sequence: {self.cornercase_seq.get_name()}")
        # Runs ONLY the chirp sequence
        await self.cornercase_seq.start(self.env.agt.sqr)
        self.logger.info(f"finished sequence: {self.cornercase_seq.get_name()}")
        self.drop_objection()   
# ---------------------------------------------------------
# Test 7: Corner Case only test (For fast debugging of corner cases and edge conditions)
# ---------------------------------------------------------                
@pyuvm.test()
class dds_fft_boundary_only_test(base_test):
    def build_phase(self):
        super().build_phase()
        self.fft_boundary_seq = dds_fft_boundary_seq.create("fft_boundary_seq")
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        await self.generate_clock()
        await self.run_initial_setup()
        self.logger.info(f"Starting sequence: {self.fft_boundary_seq.get_name()}") 
        # Runs ONLY the chirp sequence
        await self.fft_boundary_seq.start(self.env.agt.sqr)
        self.logger.info(f"finished sequence: {self.fft_boundary_seq.get_name()}")
        self.drop_objection()            
               
# ---------------------------------------------------------
# Test 8: The Full Regression (Runs everything)
# ---------------------------------------------------------                 
@pyuvm.test()
class dds_full_test(base_test):
    def build_phase(self):
        super().build_phase()
        ConfigDB.set(self, "*", "VERIF_MODE", "DDS")  # Set the verification mode to TOP for this test

        # write sequences here >>>>>>>
        self.random_seq = dds_random_seq.create("random_seq")
        self.singletone_seq = dds_singletone_seq.create("singletone_seq")
        self.cyclestress_seq = dds_cyclestress_seq.create("cyclestress_seq")
        self.chirpsweep_seq = dds_chirpsweep_seq.create("chirpsweep_seq")
        self.reset_recovery_seq = dds_reset_recovery_seq.create("dds_reset_recovery_seq")
        self.cornercase_seq= dds_cornercase_seq.create("dds_cornercase_seq")
        self.fft_boundary_seq = dds_fft_boundary_seq.create("dds_fft_boundary_seq")
        

    # run phase
    async def run_phase(self):
        test_sequences = [
            self.cornercase_seq,
            self.fft_boundary_seq,
            self.singletone_seq,
            self.cyclestress_seq,
            self.chirpsweep_seq,
            self.reset_recovery_seq,
            self.random_seq          # Brute-force random runs last!
        ]
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        # 2. Call the universal setup from base_test (Reset)
        await self.run_initial_setup()

        # 3. Run the sequences
        for seq in test_sequences:
            self.logger.info(f"Starting sequence: {seq.get_name()}")
            await seq.start(self.env.agt.sqr)
            self.logger.info(f"Finished sequence: {seq.get_name()}")

        # 4. Drop the objection to end the simulation
        self.drop_objection()
         
                

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")