"""

"""
import pyuvm
from pyuvm import *
from top_seq_item import *
import cocotb
from cocotb.triggers import RisingEdge

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
            
            # 3. Drive the stimulus into the dut_drv inputs
            self.dut_drv.rst_n.value = seq_item.rst_n
            self.dut_drv.dds_enable.value = seq_item.dds_enable
            self.dut_drv.FTW_start.value = seq_item.FTW_start
            self.dut_drv.cycles.value = seq_item.cycles
            self.dut_drv.FTW_step.value = seq_item.FTW_step

            #self.dut_drv.ofdm_in_re.value = seq_item.ofdm_in_real
            #self.dut_drv.ofdm_in_im.value = seq_item.ofdm_in_imag

            # printing the driven values for debugging
            self.logger.info(
                f"=================================================================================================================================================="
                f"Driving: rst_n={seq_item.rst_n}, dds_enable={seq_item.dds_enable}, FTW_start={seq_item.FTW_start}, cycles={seq_item.cycles}, FTW_step={seq_item.FTW_step}")
            
            # 4. Notify the sequencer that this transaction is complete
            self.seq_item_port.item_done()