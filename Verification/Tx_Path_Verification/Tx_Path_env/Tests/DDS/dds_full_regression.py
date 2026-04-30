"""
"""
import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

# Import the base test and the specific sequences
from dds_base_test import *
from dds_sequences import * 
# ---------------------------------------------------------
# Test 8: The Full Regression (Runs everything)
# ---------------------------------------------------------                 
@pyuvm.test()
class dds_full_test(dds_base_test):
    def build_phase(self):
        super().build_phase()
        ConfigDB().set(self, "*", "VERIF_MODE", "DDS")  # Set the verification mode to TOP for this test

        # write sequences here >>>>>>>
        self.random_seq = dds_random_seq.create("random_seq")
        self.singletone_seq = dds_singletone_seq.create("singletone_seq")
        self.cyclestress_seq = dds_cyclestress_seq.create("cyclestress_seq")
        self.chirpsweep_seq = dds_chirpsweep_seq.create("chirpsweep_seq")
        self.cornercase_seq= dds_cornercase_seq.create("dds_cornercase_seq")
       
        

    # run phase
    async def run_phase(self):
        test_sequences = [
            self.cornercase_seq,
            self.singletone_seq,
            self.cyclestress_seq,
            self.chirpsweep_seq,
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
            await seq.start(self.env.dds_agt.sqr)
            self.logger.info(f"Finished sequence: {seq.get_name()}")

        # 4. Drop the objection to end the simulation
        self.drop_objection()
    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")