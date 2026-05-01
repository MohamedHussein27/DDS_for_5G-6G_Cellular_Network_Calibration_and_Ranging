import cocotb
from cocotb.triggers import RisingEdge, FallingEdge
from pyuvm import *
from dds_seq_item import dds_seq_item

class dds_driver(uvm_driver): 
    def build_phase(self):
        super().build_phase()
        self.dut = ConfigDB().get(self,"","DUT")
          
    async def run_phase(self):
        # 1. Initialize hardware
        self.dut.rst_n.value = 1
        self.dut.enable.value = 0
        self.dut.FTW_start.value = 0
        self.dut.FTW_step.value = 0
        self.dut.cycles.value = 0
        
        await RisingEdge(self.dut.clk)

        # 2. Main Loop
        while True: 
            req = await self.seq_item_port.get_next_item()
            await FallingEdge(self.dut.clk) 
            
            self.dut.rst_n.value = req.rst_n
            self.logger.info(f" driver sent: rst_n={req.rst_n}, FTW_start={req.FTW_start}, FTW_step={req.FTW_step}, cycles={req.cycles}")
            if req.rst_n == 0:
                self.dut.enable.value = 0
                for _ in range(3):
                    await RisingEdge(self.dut.clk)
                self.seq_item_port.item_done()
                continue

            # Drive the configuration inputs FIRST
            self.dut.FTW_start.value = req.FTW_start
            self.dut.FTW_step.value = req.FTW_step
            self.dut.cycles.value = req.cycles
            self.dut.rst_n.value = req.rst_n
            
            # --- THE FIX: WAIT 1 CLOCK CYCLE ---
            # Allow the Verilog registers to settle before triggering the enable latch
            await RisingEdge(self.dut.clk)
            await FallingEdge(self.dut.clk)
            
            # NOW trigger the burst
            self.dut.enable.value = 1  
            
            for _ in range(req.cycles):
                await RisingEdge(self.dut.clk)
                
            await FallingEdge(self.dut.clk)
            self.dut.enable.value = 0
            
            # Pipeline Flush Delay
            for _ in range(5):
                await RisingEdge(self.dut.clk)
            
            self.seq_item_port.item_done()
        