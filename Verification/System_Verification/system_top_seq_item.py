"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: system_top_seq_item.py

    Description:
        Unified transaction-level packet for the SYSTEM_TOP block.

        This item is the currency exchanged between every component of
        the system-level UVM environment.  It merges the fields of
        top_item (TX side) and rx_item (RX side) into a single object
        so that the system monitor can stamp one transaction per clock
        and the system scoreboard can compare all three output paths
        (TX + OFDM + Radar) from the same packet.

        Field layout
        ─────────────
        ┌── System-level controls ──────────────────────────────────┐
        │  rst_n                                                     │
        │                                                            │
        ├── TX-side memory-mapped bus (from TX agent) ──────────────┤
        │  addr, wr_en, wr_data, rd_en                              │
        │  FTW_start, FTW_step (derived, not driven directly)       │
        │  dds_ready_flag, rd_data  (DUT outputs captured by mon)   │
        │  is_backdoor, backdoor_re[], backdoor_im[]                │
        │                                                            │
        ├── TX datapath outputs (monitored) ────────────────────────┤
        │  tx_valid, tx_out_real, tx_out_imag                       │
        │                                                            │
        ├── RX datapath inputs (loopback, for coverage/logging) ────┤
        │  rx_valid_in, rx_in_re, rx_in_im                         │
        │                                                            │
        ├── OFDM output path (monitored) ───────────────────────────┤
        │  ofdm_valid_out, ofdm_out_re, ofdm_out_im                │
        │                                                            │
        └── Radar output path (monitored) ──────────────────────────┘
           radar_valid_out, radar_out_re, radar_out_im

        CRV setup mirrors top_item: only bus-level signals are
        randomized; FTW registers are derived via calculate_hw_registers().
"""

import random
import math
from cocotb_coverage.crv import Randomized
from pyuvm import uvm_sequence_item


class system_top_item(uvm_sequence_item, Randomized):

    def __init__(self, name="system_top_item"):
        uvm_sequence_item.__init__(self, name)
        Randomized.__init__(self)

        # ── System Word Length constant ──────────────────────────────
        self.Fs = 491_520_000   # sampling frequency (exact integer)
        self.M  = 32            # FTW accumulator width

        # ── High-level physical parameters (used for FTW calculation) ─
        self.f0            = 0
        self.B             = 0
        self.target_cycles = 4096

        # ─────────────────────────────────────────────────────────────
        # TX-SIDE INPUTS  (driven by tx_agent driver)
        # ─────────────────────────────────────────────────────────────
        self.rst_n    = 1
        self.addr     = 0
        self.wr_en    = 0
        self.wr_data  = 0
        self.rd_en    = 0

        # Derived hardware registers (calculated, not randomized)
        self.FTW_start = 0
        self.FTW_step  = 0

        # ── Backdoor ROM write fields (tx_agent driver only) ─────────
        self.is_backdoor  = False
        self.backdoor_re  = []
        self.backdoor_im  = []

        # ─────────────────────────────────────────────────────────────
        # TX-SIDE OUTPUTS  (captured by system monitor)
        # ─────────────────────────────────────────────────────────────
        self.rd_data         = 0
        self.dds_ready_flag  = 0
        self.tx_valid        = 0
        self.tx_out_real     = 0
        self.tx_out_imag     = 0

        # ─────────────────────────────────────────────────────────────
        # RX-SIDE LOOPBACK INPUTS  (monitored, not driven externally)
        # ─────────────────────────────────────────────────────────────
        self.rx_valid_in = 0
        self.rx_in_re    = 0
        self.rx_in_im    = 0

        # ─────────────────────────────────────────────────────────────
        # OFDM OUTPUT PATH  (captured by system monitor)
        # ─────────────────────────────────────────────────────────────
        self.ofdm_valid_out = 0
        self.ofdm_out_re    = 0
        self.ofdm_out_im    = 0

        # ─────────────────────────────────────────────────────────────
        # RADAR OUTPUT PATH  (captured by system monitor)
        # ─────────────────────────────────────────────────────────────
        self.radar_valid_out = 0
        self.radar_out_re    = 0
        self.radar_out_im    = 0

        # ─────────────────────────────────────────────────────────────
        # CRV SETUP  (only bus-control signals are randomized)
        # ─────────────────────────────────────────────────────────────
        self.add_rand("rst_n",   [0, 1])
        self.add_rand("wr_en",   [0, 1])
        self.add_rand("rd_en",   [0, 1])
        self.add_rand("addr",    [0x0, 0x4, 0x8])

    # ── FTW helpers (mirrors top_item exactly) ────────────────────────
    def calculate_hw_registers(self):
        """Derive FTW_start / FTW_step from f0, B, target_cycles."""
        self.FTW_start = int(math.floor((self.f0 * (1 << self.M)) / self.Fs))
        self.FTW_step  = int(math.floor(
            (self.B * (1 << self.M)) / (self.Fs * self.target_cycles)
        ))
        mask = (1 << self.M) - 1
        self.FTW_start &= mask
        self.FTW_step  &= mask

    def calculate_chirp(self, f0, B, cycles=4096):
        """Helper to program a specific chirp without randomizing."""
        self.f0            = f0
        self.B             = B
        self.target_cycles = cycles
        self.calculate_hw_registers()

    # ── Bus helpers (mirrors top_item exactly) ────────────────────────
    def set_bus_write(self, addr):
        """Route FTW registers onto the physical bus for a write cycle."""
        self.wr_en = 1
        self.rd_en = 0
        self.addr  = addr
        if   addr == 0x0: self.wr_data = self.FTW_start
        elif addr == 0x4: self.wr_data = self.FTW_step
        elif addr == 0x8: self.wr_data = self.target_cycles
        else:             self.wr_data = 0

    def set_bus_read(self, addr=0x0):
        self.wr_en   = 0
        self.rd_en   = 1
        self.addr    = addr
        self.wr_data = 0

    def set_bus_idle(self):
        self.wr_en   = 0
        self.rd_en   = 0
        self.addr    = 0
        self.wr_data = 0

    # ── Backdoor helper ───────────────────────────────────────────────
    def set_backdoor_rom(self, re_array, im_array):
        """Pack a zero-time backdoor OFDM ROM write transaction."""
        self.is_backdoor = True
        self.backdoor_re = list(re_array)
        self.backdoor_im = list(im_array)

    # ── Randomization ─────────────────────────────────────────────────
    def randomize(self):
        """Randomize bus-level signals; derive FTW from physical limits."""
        nyquist  = self.Fs / 2
        self.f0  = random.uniform(0, nyquist / 2)
        self.B   = random.uniform(1e6, nyquist - self.f0)
        self.rst_n = random.choices([0, 1], weights=[5, 95])[0]
        self.calculate_hw_registers()

    # ── String formatting ─────────────────────────────────────────────
    def convert2string(self):
        return (
            f"{self.get_name()} "
            f"BUS[ addr={hex(self.addr)}, we={self.wr_en}, "
            f"wdata={hex(self.wr_data)}, re={self.rd_en} ] | "
            f"TX[ vld={self.tx_valid}, re={self.tx_out_real}, im={self.tx_out_imag} ] | "
            f"OFDM[ vld={self.ofdm_valid_out}, re={self.ofdm_out_re}, "
            f"im={self.ofdm_out_im} ] | "
            f"RADAR[ vld={self.radar_valid_out}, re={self.radar_out_re}, "
            f"im={self.radar_out_im} ]"
        )

    def convert2string_stimulus(self):
        return (
            f"{self.get_name()} "
            f"STIMULUS: rst_n={self.rst_n}, "
            f"addr={hex(self.addr)}, we={self.wr_en}, "
            f"wdata={hex(self.wr_data)}, re={self.rd_en} "
            f"[FTW_start={hex(self.FTW_start)}, FTW_step={hex(self.FTW_step)}]"
        )

    def convert2string_output(self):
        return (
            f"{self.get_name()} "
            f"OUTPUT: "
            f"tx_vld={self.tx_valid}, tx_re={self.tx_out_real}, tx_im={self.tx_out_imag} | "
            f"ofdm_vld={self.ofdm_valid_out}, ofdm_re={self.ofdm_out_re}, "
            f"ofdm_im={self.ofdm_out_im} | "
            f"radar_vld={self.radar_valid_out}, radar_re={self.radar_out_re}, "
            f"radar_im={self.radar_out_im}"
        )
