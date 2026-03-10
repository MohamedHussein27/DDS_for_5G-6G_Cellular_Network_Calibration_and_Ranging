from cocotb.triggers import *
import cocotb

from transaction import Transaction

class Driver:

    def __init__(self):
        self.gen2drv = None
        self.drv_rqt = Event()
        self.tr_drv = Transaction()

    async def run(self, dut):
        
        while True:

            self.tr_drv = await self.gen2drv.get()
            await FallingEdge(dut.clk)
            
            cocotb.log.info("[Driver] begins.")

            dut.rst_n.value = self.tr_drv.rst_n
            dut.a.value     = self.tr_drv.a
            dut.b.value     = self.tr_drv.b
            dut.op.value    = self.tr_drv.op

            
            self.drv_rqt.set()
            """if self.drv_rqt is not None:
                self.drv_rqt.set()
                self.drv_rqt.clear()"""