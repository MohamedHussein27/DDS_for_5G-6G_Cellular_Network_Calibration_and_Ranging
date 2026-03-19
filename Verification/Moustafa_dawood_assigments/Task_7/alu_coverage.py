import pyuvm
from pyuvm import *

class alu_coverage(uvm_component):
    def build_phase(self):
        super().build_phase()
        # 1. TRADITION: Creation Phase
        self.cov_fifo = uvm_tlm_analysis_fifo("cov_fifo", self)

    def connect_phase(self):
        super().connect_phase()
        # 2. TRADITION: Wiring Phase
        self.analysis_export = self.cov_fifo.analysis_export

    async def run_phase(self):
        while True:
            item = await self.cov_fifo.get()
            pass