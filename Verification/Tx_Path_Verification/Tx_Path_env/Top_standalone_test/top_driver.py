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
            
            # 2. Synchronize with the clock edge before driving the signals
            await FallingEdge(self.dut_drv.clk)
            
            # 3. Drive the stimulus into the memory-mapped bus inputs
            self.dut_drv.rst_n.value   = seq_item.rst_n
            self.dut_drv.addr.value    = seq_item.addr
            self.dut_drv.wr_en.value   = seq_item.wr_en
            self.dut_drv.wr_data.value = seq_item.wr_data
            self.dut_drv.rd_en.value   = seq_item.rd_en

            # printing the driven values for debugging
            self.logger.info(
                f"====================================================================================================\n"
                f"Driving BUS: rst_n={seq_item.rst_n}, addr={hex(seq_item.addr)}, we={seq_item.wr_en}, "
                f"wdata={hex(seq_item.wr_data)}, re={seq_item.rd_en}"
            )
            
            # 4. Notify the sequencer that this transaction is complete
            self.seq_item_port.item_done()