import pyuvm
from pyuvm import *
import cocotb
from cocotb.triggers import FallingEdge

class alu_driver(uvm_driver):
    def build_phase(self):
        super().build_phase()

    async def run_phase(self):
        while True:
            item = await self.seq_item_port.get_next_item()
            
            # Driving actual Verilog pins
            self.dut.op.value = item.op
            self.dut.a.value = item.a
            self.dut.b.value = item.b
            self.dut.rst_n.value = item.rst_n
            
            await FallingEdge(self.dut.clk)
            
            self.seq_item_port.item_done()