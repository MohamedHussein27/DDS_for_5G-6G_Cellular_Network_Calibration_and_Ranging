"""

"""
import pyuvm
from pyuvm import *
from seq_item import *
import cocotb
from cocotb.triggers import RisingEdge

class top_driver(uvm_driver):
    def build_phase(self):
        # Fetch the DUT from the UVM Configuration Database
        self.dut = ConfigDB().get(self, "", "DUT")
        
    async def run_phase(self):
        while True:
            # 1. Get the next transaction item from the sequencer
            seq_item = await self.seq_item_port.get_next_item()
            
            # 2. Synchronize with the clock edge before driving the signals
            await FallingEdge(self.dut.clk)
            
            # 3. Drive the stimulus into the DUT inputs
            self.dut.rst_n.value = seq_item.rst_n
            self.dut.dds_enable.value = seq_item.dds_enable
            self.dut.FTW_start.value = seq_item.FTW_start
            self.dut.cycles.value = seq_item.cycles
            self.dut.FTW_step.value = seq_item.FTW_step
            self.dut.ofdm_rd_en.value = seq_item.ofdm_rd_en
            self.dut.ofdm_in_real.value = seq_item.ofdm_in_real
            self.dut.ofdm_in_imag.value = seq_item.ofdm_in_imag
            
            # 4. Notify the sequencer that this transaction is complete
            self.seq_item_port.item_done()