import cocotb
from cocotb.triggers import RisingEdge
from transaction import ALUTransaction
from cocotb.triggers import Timer
class ALUMonitor:

    def __init__(self, dut, sb, cov):
        self.dut = dut
        self.sb = sb
        self.cov = cov

    async def run(self):

        while True:

            await RisingEdge(self.dut.clk)
            
            await Timer(1, units="ns")
            
            item = ALUTransaction(
            self.dut.a.value.to_unsigned(),
            self.dut.b.value.to_unsigned(),
            self.dut.op.value.to_unsigned()
)
            
            
            
            c_val = int(self.dut.c.value)
            out_val = self.dut.out.value.to_unsigned()
            
            result = (c_val << 4) | out_val
            item.result = result

            self.sb.check(item)
            self.cov.sample(item)