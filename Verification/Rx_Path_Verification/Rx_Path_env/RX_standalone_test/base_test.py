"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: base_test.py

    Description:
        The base test serves as the foundational pyuvm test class for the TX wrapper
        verification environment. It establishes the testbench hierarchy, instantiates
        the environment components, generates the 491.52 MHz system clock, and defines
        the standard execution flow (reset only and creating the clock) that all
        extended tests will inherit.
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm

from environment import *

class base_test(uvm_test):

    # FIX A: was `_init_` (single underscores) — Python never called it.
    def __init__(self, name, parent):
        super().__init__(name, parent)

    def build_phase(self):

        # Environment
        self.env = Environment.create("env", self)

        # dut handle — the top-level RX_TOP cocotb object
        self.dut = cocotb.top

        # ── DUT handle distribution ──────────────────────────────────────
        # top_agt and top_sb/sub get the full RX_TOP handle (for rx_in, etc.)
        ConfigDB().set(self, "*", "DUT", self.dut)

        # fft_agt monitor reads u_rx_fft submodule ports (valid_in, out_real …)
        # Those ARE input ports driven externally by rx_valid_in/rx_in_re/re_in_im,
        # so cocotb can read them correctly from the submodule handle.
        ConfigDB().set(self, "env.fft_agt.*", "DUT", self.dut.u_rx_fft)

        # FIX B: ifft_agt monitor must get the TOP-level handle, NOT u_rx_ifft.
        #
        # WHY: u_rx_ifft.valid_in is driven by the internal wire `mult_valid`
        # inside RX_TOP.  Cocotb cannot back-read a submodule INPUT port through
        # the hierarchy boundary — it always returns 0, so the scoreboard sees
        # 0/0 samples even though the RTL is working.
        #
        # The ifft_monitor has been updated to read the driving wires directly
        # from the top level:
        #   valid_in  <- dut.mult_valid
        #   in_real   <- dut.mult_re >> 3   (replicating the RTL >>>3 shift)
        #   in_imag   <- dut.mult_im >> 3
        #   valid_out <- dut.ifft_valid
        #   out_real  <- dut.ifft_re
        #   out_imag  <- dut.ifft_im
        ConfigDB().set(self, "env.ifft_agt.*", "DUT", self.dut)   # <-- top-level handle

        # ── Agent activity ───────────────────────────────────────────────
        ConfigDB().set(self, "env.top_agt",  "is_active", uvm_active_passive_enum.UVM_ACTIVE)
        ConfigDB().set(self, "env.fft_agt",  "is_active", uvm_active_passive_enum.UVM_PASSIVE)
        ConfigDB().set(self, "env.ifft_agt", "is_active", uvm_active_passive_enum.UVM_PASSIVE)

    async def generate_clock(self):
        self.clk = Clock(self.dut.clk, 2034, units="ps")  # 491.52 MHz
        await cocotb.start(self.clk.start())

    def final_phase(self):
        self.logger.info(f"**** End of {self.get_type_name()} ****")