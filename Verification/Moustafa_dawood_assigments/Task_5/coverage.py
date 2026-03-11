from cocotb_coverage.coverage import CoverPoint, CoverCross


class ALUCoverage:

    def __init__(self):

        self.mon2cov = None

    @CoverPoint("alu.a",
                xf=lambda self, item: item.a,
                bins=list(range(16)))

    @CoverPoint("alu.b",
                xf=lambda self, item: item.b,
                bins=list(range(16)))

    @CoverPoint("alu.op",
                xf=lambda self, item: item.op,
                bins=[0,1,2,3])

    @CoverCross("alu.cross",
                items=["alu.a","alu.b","alu.op"])

    def sample(self, item):
        pass

    async def run(self):

        while True:

            item = await self.mon2cov.get()

            self.sample(item)