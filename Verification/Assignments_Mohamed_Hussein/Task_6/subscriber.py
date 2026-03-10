from cocotb_coverage.coverage import CoverPoint, CoverCross, coverage_section
from cocotb.queue import *
from cocotb_coverage.coverage import coverage_db
import cocotb

class Subscriber:

    def __init__(self):

        self.mon2cov = Queue()
        self.cov_item = None


    # -----------------------------
    # COVERAGE DEFINITIONS
    # -----------------------------

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


    # -----------------------------
    # RUN
    # -----------------------------

    async def run(self):

        while True:

            item = await self.mon2cov.get()

            self.cov_item = item

            # sample coverage
            self.sample(item)


    # -----------------------------
    # REPORT
    # -----------------------------

    def report(self):
        

        cocotb.log.info("========= FUNCTIONAL COVERAGE REPORT =========")
        # Pass the cocotb logger so it prints perfectly in the transcript
        coverage_db.report_coverage(cocotb.log.info, bins=True)
        cocotb.log.info("==============================================")