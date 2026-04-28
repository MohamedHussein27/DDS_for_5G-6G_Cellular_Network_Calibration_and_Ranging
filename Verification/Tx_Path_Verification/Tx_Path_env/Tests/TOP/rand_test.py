"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: rand_test.py

    Description:
        This test verifies the TX datapath under extreme, unconstrained randomized 
        conditions. It dynamically selects and executes any available sequence at 
        random to stress the datapath, pipeline routing, and multiplexer toggling. 
        It ensures that unpredictable switching between null frames, resets, single-band 
        operation, and edge-case guardband constraints (strictly 200KHz to 210KHz) 
        does not corrupt the Q5.26 fixed-point precision or cause pipeline deadlocks.
"""
import cocotb 
from cocotb.triggers import * 
from cocotb.clock import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging
import random

# Import the base test and ALL specific sequences across the 5 groups
from base_test import base_test
from null_sequences import null_frame_seq, ofdm_only_seq, chirp_only_seq
from single_band_sequences import single_chirp_tone_seq, full_dds_sweep_seq
from word_format_sequences import word_format_normal_seq, overflow_saturation_seq, valid_out_gating_seq
from reset_sequences import reset_before_frame_seq, reset_mid_frame_seq, multiple_resets_seq
from regression_sequences import back_to_back_frames_seq, long_gap_frames_seq, random_frames_regression_seq

@pyuvm.test()
class rand_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Create a comprehensive pool of ALL 14 diverse sequences to randomly choose from
        self.seq_pool = [
            # Group A: Null & Boundary
            null_frame_seq.create("seq_null_both"),
            ofdm_only_seq.create("seq_ofdm_only"),
            chirp_only_seq.create("seq_chirp_only"),
            
            # Group B: Single-Band Operation
            single_chirp_tone_seq.create("seq_single_tone"),
            full_dds_sweep_seq.create("seq_full_sweep"),
            
            # Group C: Word Format & Overflow
            word_format_normal_seq.create("seq_word_format"),
            overflow_saturation_seq.create("seq_overflow"),
            valid_out_gating_seq.create("seq_valid_gating"),

            # Group F: Reset & Control
            reset_before_frame_seq.create("seq_reset_before"),
            reset_mid_frame_seq.create("seq_reset_mid"),
            multiple_resets_seq.create("seq_reset_multi"),

            # Group H: Regression & Stress
            back_to_back_frames_seq.create("seq_back_to_back"),
            long_gap_frames_seq.create("seq_long_gap"),
            random_frames_regression_seq.create("seq_full_regress")
        ]

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        # 2. Call the universal setup from base_test (Reset)
        await self.run_initial_setup()

        # 3. Random Sequence Execution Loop
        num_iterations = 50  # Define how many random sequences to fire back-to-back
        self.logger.info(f"--- Firing {num_iterations} completely random sequences ---")
        
        for i in range(num_iterations):
            # Randomly select a sequence from the entire verification suite pool
            selected_seq = random.choice(self.seq_pool)
            self.logger.info(f"[Iteration {i+1}/{num_iterations}] Randomly selected: {selected_seq.get_name()}")
            
            # Execute the randomly selected sequence
            await selected_seq.start(self.env.agt.sqr)

        # 4. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Extreme dynamic sequence switching across all 14 test cases")
        self.logger.info("---------------------------------------------------------")