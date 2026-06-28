import pyuvm
from pyuvm import *
from top_seq_item import *
import cocotb
from cocotb.triggers import *

class top_driver(uvm_driver):
    def build_phase(self):
        # Fetch the DUT from the UVM Configuration Database
        self.dut_drv = ConfigDB().get(self, "", "DUT")
        
    async def run_phase(self):
        while True:
            # 1. Get the next transaction item from the sequencer
            seq_item = await self.seq_item_port.get_next_item()
            
            # ==========================================================
            # INTERCEPT BACKDOOR TRANSACTIONS (Zero-Time Execution)
            # ==========================================================
            if seq_item.is_backdoor:

                rom_re = self.dut_drv.u_ofdm_rom.rom_real
                rom_im = self.dut_drv.u_ofdm_rom.rom_imag
                
                # Write the payload arrays into the simulation memory
                for idx in range(len(seq_item.backdoor_re)):
                    rom_re[idx].value = seq_item.backdoor_re[idx]
                    rom_im[idx].value = seq_item.backdoor_im[idx]
                    
                self.logger.info(f"Driver executed ZERO-TIME backdoor ROM write ({len(seq_item.backdoor_re)} entries).")
                self.seq_item_port.item_done()
                continue # Skip the rest of the loop! Do not wiggle pins!
            
            # ==========================================================
            # STANDARD PHYSICAL BUS WIGGLING
            # ==========================================================
            await FallingEdge(self.dut_drv.clk)
            
            self.dut_drv.rst_n.value   = seq_item.rst_n
            self.dut_drv.addr.value    = seq_item.addr
            self.dut_drv.wr_en.value   = seq_item.wr_en
            self.dut_drv.wr_data.value = seq_item.wr_data
            self.dut_drv.rd_en.value   = seq_item.rd_en

            # printing the driven values for debugging
            self.logger.debug(
                f"====================================================================================================\n"
                f"Driving BUS: rst_n={seq_item.rst_n}, addr={hex(seq_item.addr)}, we={seq_item.wr_en}, "
                f"wdata={hex(seq_item.wr_data)}, re={seq_item.rd_en}"
            )
            
            # 4. Notify the sequencer that this transaction is complete
            self.seq_item_port.item_done()