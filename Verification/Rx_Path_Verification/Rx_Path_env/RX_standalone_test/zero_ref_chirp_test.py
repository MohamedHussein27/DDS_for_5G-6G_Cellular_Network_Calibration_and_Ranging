"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: zero_ref_chirp_test.py

    Description:
        This test verifies the RX datapath under TC-RX-002:
        Zero Reference, Chirp Signal In.

        The reference RAM is never written (stays at power-on zero).
        A chirp signal with energy placed only in radar frequency bins
        (bins 2048-4095 after the mid-path bit-reversal) is driven
        into rx_in for one full 4096-sample frame.

    STIMULUS:
        Apply chirp-only signal (radar bins only).
        Reference RAM = zeros (never written).
        Capture radar_out and ofdm_out.

    EXPECTED / PASS:
        radar_out = 0 on all 2048 output samples.
            (conj_multiply × zero_reference → zero regardless of input)
        ofdm_out matches the golden FFT model.
            (OFDM bins 0-2047 contain zero because the chirp was placed
             only in the radar half of the spectrum)

    Feature: F-RX-02
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from pyuvm import *
import pyuvm
import logging

from base_test                    import base_test
from top_seq_item                  import rx_item
from zero_ref_chirp_seq import *


@pyuvm.test()
class zero_ref_chirp_test(base_test):
    """
    TC-RX-002: Zero Reference, Chirp Signal In.
    Verifies that the radar output is exactly zero when the reference RAM
    is never populated, and that the OFDM output matches the golden model.
    """

    def build_phase(self):
        super().build_phase()

        # Override agent activity for RX TOP:
        #   top_agt  → ACTIVE  (drives rx_valid_in + ref_wr_en)
        #   fft_agt  → PASSIVE (monitors FFT boundary)
        #   ifft_agt → PASSIVE (monitors IFFT boundary)
        ConfigDB().set(self, "env.top_agt",  "is_active",
                       uvm_active_passive_enum.UVM_ACTIVE)
        ConfigDB().set(self, "env.fft_agt",  "is_active",
                       uvm_active_passive_enum.UVM_PASSIVE)
        ConfigDB().set(self, "env.ifft_agt", "is_active",
                       uvm_active_passive_enum.UVM_PASSIVE)

        # Set full-system verification mode
        ConfigDB().set(self, "env", "VERIF_MODE", "TOP")

        # Instantiate the TC-RX-002 sequence
        self.seq_tc_rx_002 = zero_ref_chirp_seq.create(
            "seq_tc_rx_002"
        )

        # Configure: run one frame (extend to stress-test by setting > 1)
        self.seq_tc_rx_002.num_frames = 1

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )

        # 1. Start the 491.52 MHz system clock
        await self.generate_clock()

        # 2. Execute TC-RX-002
        self.logger.info(
            "--- Executing TC-RX-002: Zero Reference, Chirp Signal In ---"
        )
        await self.seq_tc_rx_002.start(self.env.top_agt.sqr)

        self.drop_objection()

    def report_phase(self):
        self.logger.info(
            "---------------------------------------------------------"
        )
        self.logger.info(f" [TEST REPORT]   : {self.get_type_name()}")
        self.logger.info(
            " [TC ID]         : TC-RX-002"
        )
        self.logger.info(
            " [STIMULUS]      : Chirp signal (radar bins only). "
            "Reference RAM = zeros (never written)."
        )
        self.logger.info(
            " [PASS CRITERIA] : "
            "radar_out = 0 on all 2048 samples. "
            "ofdm_out matches golden FFT model (zero OFDM bins)."
        )
        self.logger.info(
            " [FEATURE]       : F-RX-02 — "
            "Radar path: conjugate-multiply drives output to zero "
            "when reference is unpopulated."
        )
        self.logger.info(
            "---------------------------------------------------------"
        )