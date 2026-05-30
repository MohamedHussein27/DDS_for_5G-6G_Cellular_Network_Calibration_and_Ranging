"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: system_top_monitor.py

    Description:
        Monitor for the SYSTEM_TOP UVM agent.

        Observes ALL DUT boundary signals on every rising clock edge and
        packs them into a single system_top_item transaction that is
        broadcast through mon_ap to any connected subscriber or
        scoreboard.

        Captured signal groups
        ──────────────────────
        ① TX-side bus inputs
              rst_n, addr, wr_en, wr_data, rd_en
        ② TX-side bus outputs
              rd_data, dds_ready_flag
        ③ TX datapath outputs
              tx_valid, tx_out_re → tx_out_real, tx_out_im → tx_out_imag
        ④ RX loopback inputs (internally wired in RTL, monitored for coverage)
              rx_valid_in (= rf_valid from TX), rx_in_re, rx_in_im
        ⑤ OFDM RX output
              ofdm_valid_out, ofdm_out_re, ofdm_out_im
        ⑥ Radar RX output
              radar_valid_out, radar_out_re, radar_out_im

        Unsafe RTL states ('X', 'Z') are handled defensively: 1-bit
        signals default to 0, multi-bit signed signals default to 0.

        The monitor does NOT inspect backdoor transactions — those are
        purely software constructs that never appear on RTL wires.
"""

import pyuvm
from pyuvm import *
import cocotb
from cocotb.triggers import RisingEdge, ReadOnly

from system_top_seq_item import system_top_item


class system_top_monitor(uvm_monitor):

    def build_phase(self):
        self.mon_ap  = uvm_analysis_port("mon_ap", self)
        self.dut_mon = ConfigDB().get(self, "", "DUT")

    async def run_phase(self):

        # ── Safe-read helpers ────────────────────────────────────────
        def read_1bit(sig):
            try:
                return int(sig.value)
            except ValueError:
                return 0

        def read_signed(sig):
            try:
                return sig.value.signed_integer
            except ValueError:
                return 0

        def read_unsigned(sig):
            try:
                return int(sig.value)
            except ValueError:
                return 0

        # Skip the very first rising edge (signals may be undefined)
        await RisingEdge(self.dut_mon.clk)

        while True:
            await RisingEdge(self.dut_mon.clk)
            await ReadOnly()

            item = system_top_item("mon_item")

            # ── ① TX bus inputs ─────────────────────────────────────
            item.rst_n    = read_1bit    (self.dut_mon.rst_n)
            item.addr     = read_unsigned(self.dut_mon.addr)
            item.wr_en    = read_1bit    (self.dut_mon.wr_en)
            item.wr_data  = read_unsigned(self.dut_mon.wr_data)
            item.rd_en    = read_1bit    (self.dut_mon.rd_en)

            # ── ② TX bus outputs ─────────────────────────────────────
            item.rd_data        = read_unsigned(self.dut_mon.rd_data)
            item.dds_ready_flag = read_1bit    (self.dut_mon.dds_ready_flag)

            # ── ③ TX datapath outputs ────────────────────────────────
            item.tx_valid    = read_1bit    (self.dut_mon.tx_valid)
            item.tx_out_real = read_signed  (self.dut_mon.tx_out_re)
            item.tx_out_imag = read_signed  (self.dut_mon.tx_out_im)

            # ── ④ RX loopback inputs (internal wires observed passively)
            # These wires exist inside the DUT; we read them for coverage.
            # If the hierarchy is not accessible, default gracefully to 0.
            try:
                item.rx_valid_in = read_1bit (self.dut_mon.u_rx.rx_valid_in)
                item.rx_in_re    = read_signed(self.dut_mon.u_rx.rx_in_re)
                item.rx_in_im    = read_signed(self.dut_mon.u_rx.rx_in_im)
            except AttributeError:
                item.rx_valid_in = 0
                item.rx_in_re    = 0
                item.rx_in_im    = 0

            # ── ⑤ OFDM RX output ─────────────────────────────────────
            item.ofdm_valid_out = read_1bit  (self.dut_mon.ofdm_valid_out)
            item.ofdm_out_re    = read_signed(self.dut_mon.ofdm_out_re)
            item.ofdm_out_im    = read_signed(self.dut_mon.ofdm_out_im)

            # ── ⑥ Radar RX output ─────────────────────────────────────
            item.radar_valid_out = read_1bit  (self.dut_mon.radar_valid_out)
            item.radar_out_re    = read_signed(self.dut_mon.radar_out_re)
            item.radar_out_im    = read_signed(self.dut_mon.radar_out_im)

            # ── Log and broadcast ────────────────────────────────────
            self.logger.info(
                f"SYS Monitor: "
                f"rst={item.rst_n}, addr={hex(item.addr)}, "
                f"we={item.wr_en}, wdata={hex(item.wr_data)}, "
                f"re={item.rd_en} | "
                f"rdata={hex(item.rd_data)}, rdy={item.dds_ready_flag} | "
                f"tx_vld={item.tx_valid} | "
                f"ofdm_vld={item.ofdm_valid_out} | "
                f"radar_vld={item.radar_valid_out}"
            )

            self.mon_ap.write(item)
