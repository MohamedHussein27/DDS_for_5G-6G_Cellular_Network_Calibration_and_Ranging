import cocotb
from cocotb.triggers import RisingEdge, FallingEdge
from pyuvm import *
from dds_seq_item import dds_seq_item

class dds_driver(uvm_driver): 
    def build_phase(self):
        super().build_phase()
        self.dut_drv = ConfigDB().get(self,"","DDS_DUT")
          
    async def run_phase(self):
        # 1. Initialize hardware
        self.dut_drv.rst_n.value = 1
        self.dut_drv.enable.value = 0
        self.dut_drv.FTW_start.value = 0
        self.dut_drv.FTW_step.value = 0
        self.dut_drv.cycles.value = 0
        
        await RisingEdge(self.dut_drv.clk)

        # 2. Main Loop
        while True: 
            req = await self.seq_item_port.get_next_item()
            await FallingEdge(self.dut_drv.clk) 
            
            self.dut_drv.rst_n.value = req.rst_n
            
            if req.rst_n == 0:
                self.dut_drv.enable.value = 0
                for _ in range(3):
                    await RisingEdge(self.dut_drv.clk)
                self.seq_item_port.item_done()
                continue

            # Drive the configuration inputs FIRST
            self.dut_drv.FTW_start.value = req.FTW_start
            self.dut_drv.FTW_step.value = req.FTW_step
            self.dut_drv.cycles.value = req.cycles
            self.dut_drv.rst_n.value = req.rst_n
            
            # --- THE FIX: WAIT 1 CLOCK CYCLE ---
            # Allow the Verilog registers to settle before triggering the enable latch
            await RisingEdge(self.dut_drv.clk)
            await FallingEdge(self.dut_drv.clk)
            
            # NOW trigger the burst
            self.dut_drv.enable.value = 1  
            
            for _ in range(req.cycles):
                await RisingEdge(self.dut_drv.clk)
                
            await FallingEdge(self.dut_drv.clk)
            self.dut_drv.enable.value = 0
            
            # Pipeline Flush Delay
            for _ in range(5):
                await RisingEdge(self.dut_drv.clk)
            
            self.seq_item_port.item_done()