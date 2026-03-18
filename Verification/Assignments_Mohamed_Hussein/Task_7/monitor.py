import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

from sequence_item import *

class Monitor(uvm_monitor):

    def __init__(self, name, parent):
        super().__init__(name, parent)

    # building share point
    def build_phase(self):
        self.mon_ap = uvm_analysis_port.create("mon_ap", self)
        # getting the dut handle and assign it to the monitor dut
        self.dut_mon= ConfigDB().get(self,"","DUT")
    

    async def run_phase(self):
        await RisingEdge(self.dut_mon.clk)
        while True:
            rsp_seq_item = Sequence_Item.create("rsp_seq_item")

            await RisingEdge(self.dut_mon.clk)
            await ReadOnly()
            # driving
            rsp_seq_item.rst_n = self.dut_mon.rst_n.value
            rsp_seq_item.a     = self.dut_mon.a.value
            rsp_seq_item.b     = self.dut_mon.b.value
            rsp_seq_item.op    = self.dut_mon.op.value
            rsp_seq_item.out   = self.dut_mon.out.value
            rsp_seq_item.c     = self.dut_mon.c.value

            """self.logger.info(f"[MOnitor] receiving: rst_n={rsp_seq_item.rst_n}, "
                         f"a={rsp_seq_item.a}, b={rsp_seq_item.b}, op={rsp_seq_item.op},"
                         f"out={rsp_seq_item.out}, c={rsp_seq_item.c}")"""

            # broadcasting the data to be collected by the sub and sb
            self.mon_ap.write(rsp_seq_item)
  