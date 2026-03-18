import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from cocotb_coverage.coverage import CoverPoint, CoverCross, coverage_section
from cocotb.queue import *
from cocotb_coverage.coverage import coverage_db
from pyuvm import * 
import pyuvm
import logging

from sequence_item  import * 

class Subscriber(uvm_subscriber):

    def __init__(self, name, parent):
        super().__init__(name, parent)

    def build_phase(self):
        self.sub_export = uvm_analysis_export("sub_export", self)
        self.sub_fifo = uvm_tlm_analysis_fifo("sub_fifo", self)

    # must be implemented
    def write(self, item):
        self.logger.debug(f"Subscriber received")

    @CoverPoint(
        "alu_cov.A_CP",
        xf=lambda tr: tr.a,
        bins=[0, (8,14), 15]
    )
    @CoverPoint(
        "alu_cov.B_CP",
        xf=lambda tr: tr.b,
        bins=[0, (8,14), 15]
    )
    @CoverPoint(
        "alu_cov.OP_CP",
        xf=lambda tr: tr.op,
        bins=[0,1,2,3]
    )
    @CoverPoint(
        "alu_cov.OUT_CP",
        xf=lambda tr: tr.out
    )
    @CoverPoint(
        "alu_cov.C_CP",
        xf=lambda tr: tr.c
    )

    @CoverCross(
        "alu_cov.BOUNDARY_ADD_C",
        items=["alu_cov.A_CP", "alu_cov.B_CP", "alu_cov.OP_CP"],
        ign_bins=[]
    )

    @CoverCross(
        "alu_cov.CARRY_C",
        items=["alu_cov.A_CP", "alu_cov.B_CP", "alu_cov.OP_CP"],
        ign_bins=[]
    )

    def sample(self, tr):
        pass

    async def run_phase(self):
        while True:
            seq_item_sub = await self.sub_fifo.get()
            self.sample(seq_item_sub)