"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: system_top_driver.py

    Description:
        Driver for the SYSTEM_TOP UVM agent.

        Receives system_top_item transactions from the sequencer and
        drives the physical pins of the SYSTEM_TOP DUT.

        The SYSTEM_TOP has the TX memory-mapped bus as its only
        *externally driven* input — the RX datapath is internally wired
        to the TX output, so the driver never touches RX pins directly.

        Two transaction modes are supported, exactly mirroring
        tx_top_driver.py:

        1. BACKDOOR (is_backdoor = True)
           Zero-simulation-time write to the internal OFDM ROM through
           the cocotb hierarchical handle.  Control is returned to the
           sequencer immediately with no clock edge consumed.

        2. NORMAL (is_backdoor = False)
           Physical bus wiggling on the falling edge of clk (setup/hold
           safe for the DUT's posedge logic).  Drives:
               rst_n, addr, wr_en, wr_data, rd_en
"""

import pyuvm
from pyuvm import *
import cocotb
from cocotb.triggers import FallingEdge

from system_top_seq_item import system_top_item


class system_top_driver(uvm_driver):

    def build_phase(self):
        # Fetch the DUT handle from the UVM Configuration Database.
        # The testbench registers it with key "DUT" at the top level.
        self.dut_drv = ConfigDB().get(self, "", "DUT")

    async def run_phase(self):
        while True:
            # ── 1. Fetch the next transaction ─────────────────────────
            seq_item = await self.seq_item_port.get_next_item()

            # ── 2. BACKDOOR path: zero-time ROM write ─────────────────
            if seq_item.is_backdoor:
                rom_re = self.dut_drv.u_tx.u_ofdm_rom.rom_real
                rom_im = self.dut_drv.u_tx.u_ofdm_rom.rom_imag

                for idx in range(len(seq_item.backdoor_re)):
                    rom_re[idx].value = seq_item.backdoor_re[idx]
                    rom_im[idx].value = seq_item.backdoor_im[idx]

                self.logger.info(
                    f"SYS Driver: BACKDOOR ROM write "
                    f"({len(seq_item.backdoor_re)} entries) — zero time."
                )
                self.seq_item_port.item_done()
                continue   # skip normal bus wiggling

            # ── 3. NORMAL path: physical bus wiggling ─────────────────
            # Drive on the falling edge so setup/hold times are met
            # for the DUT's rising-edge-triggered logic.
            await FallingEdge(self.dut_drv.clk)

            self.dut_drv.rst_n.value   = seq_item.rst_n
            self.dut_drv.addr.value    = seq_item.addr
            self.dut_drv.wr_en.value   = seq_item.wr_en
            self.dut_drv.wr_data.value = seq_item.wr_data
            self.dut_drv.rd_en.value   = seq_item.rd_en

            self.logger.info(
                f"SYS Driver: "
                f"rst_n={seq_item.rst_n}, "
                f"addr={hex(seq_item.addr)}, "
                f"we={seq_item.wr_en}, "
                f"wdata={hex(seq_item.wr_data)}, "
                f"re={seq_item.rd_en}"
            )

            # ── 4. Notify sequencer: transaction complete ──────────────
            self.seq_item_port.item_done()
