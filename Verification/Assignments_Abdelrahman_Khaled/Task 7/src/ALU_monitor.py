from pyuvm import *
from cocotb.triggers import RisingEdge
from ALU_seq_item import ALU_seq_item


class ALU_monitor(uvm_monitor):
    def build_phase(self):
        super().build_phase()
        self.monitor_ap = uvm_analysis_port("monitor_ap", self)
        self.cfg = ConfigDB().get(self, "", "CFG")
        self.vif = self.cfg.vif

    async def run_phase(self):
        while True:
            await RisingEdge(self.vif.clk)
            item_monitor = ALU_seq_item("item_monitor")

            # Try to read all signals and convert them to integers
            try:
                item_monitor.a = self.vif.a.value.integer
                item_monitor.b = self.vif.b.value.integer
                item_monitor.op = self.vif.op.value.integer
                item_monitor.out = self.vif.out.value.integer
                item_monitor.c = self.vif.c.value.integer

            # If cocotb hits an 'X' or 'Z', it throws a ValueError.
            # We catch it and set safe default values for that clock cycle!
            except ValueError:
                item_monitor.a = 0
                item_monitor.b = 0
                item_monitor.op = 0
                item_monitor.out = 0
                item_monitor.c = 0

            self.monitor_ap.write(item_monitor)
