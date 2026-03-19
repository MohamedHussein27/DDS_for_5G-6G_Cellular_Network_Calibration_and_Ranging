import pyuvm
from pyuvm import *
import cocotb
from cocotb.triggers import RisingEdge, Timer
from alu_sequence_item import alu_sequence_item
from alu_pkg import opcode_t  # <-- Added to bring the Enum into the Monitor!

class alu_monitor(uvm_monitor):
    def build_phase(self):
        super().build_phase()
        self.mon_ap = uvm_analysis_port("mon_ap", self)

    async def run_phase(self):
        while True:
            await RisingEdge(self.dut.clk)
            await Timer(1, unit="ns") 
            
            item = alu_sequence_item("item")
            
            # Reading the ACTUAL Verilog ports
            try:
                # Cast the raw integer back into our Enum type so .name works!
                raw_op = int(self.dut.op.value)
                item.op = opcode_t(raw_op)
            except ValueError:
                
                item.op = opcode_t.ADD_op 
                
            item.a = int(self.dut.a.value)
            item.b = int(self.dut.b.value)
            item.rst_n = int(self.dut.rst_n.value)
            
            # Reconstruct the 5-bit result from 'c' (carry) and 'out'
            c_val = int(self.dut.c.value)
            out_val = int(self.dut.out.value)
            item.result = (c_val << 4) | out_val
            
            self.logger.debug(f"Monitor sampled: {item.convert2string()}")
            self.mon_ap.write(item)