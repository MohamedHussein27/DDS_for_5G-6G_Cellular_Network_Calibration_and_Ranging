from pyuvm import *
from cocotb.triggers import FallingEdge


class ALU_driver(uvm_driver):
    def build_phase(self):
        super().build_phase()
        self.cfg = ConfigDB().get(self, "", "CFG")
        self.vif = self.cfg.vif

    async def run_phase(self):
        # --- Drive default values at time 0 to clear the blue 'Z' lines ---
        self.vif.a.value = 0
        self.vif.b.value = 0
        self.vif.op.value = 0

        while True:
            # Wait for the next transaction from the sequence
            item_dr = await self.seq_item_port.get_next_item()

            # Drive the transaction values
            self.vif.a.value = item_dr.a
            self.vif.b.value = item_dr.b
            self.vif.op.value = item_dr.op

            # Wait for the clock edge, then tell the sequencer we are done
            await FallingEdge(self.vif.clk)
            self.seq_item_port.item_done()
            self.logger.debug(f"Driven item: {item_dr}")
