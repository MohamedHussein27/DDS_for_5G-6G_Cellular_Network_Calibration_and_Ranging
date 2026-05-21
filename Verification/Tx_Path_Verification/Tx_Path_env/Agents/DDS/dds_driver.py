import cocotb
from cocotb.triggers import RisingEdge, FallingEdge
from pyuvm import *
from dds_seq_item import dds_seq_item

class dds_driver(uvm_driver): 
    def build_phase(self):
        super().build_phase()
        # Retrieve the DUT handle from the UVM configuration database
        self.dut = ConfigDB().get(self,"","DUT")
          
    async def run_phase(self):
        # 1. Initialize hardware: Set default/idle states for all inputs
        self.dut.rst_n.value = 1
        self.dut.enable.value = 0
        self.dut.FTW_start.value = 0
        self.dut.FTW_step.value = 0
        self.dut.cycles.value = 0
        
        # Wait for the first clock edge before starting the main driver loop
        await RisingEdge(self.dut.clk)

        # 2. Main Loop: Fetch and drive sequence items continuously
        while True: 
            # Get the next transaction item from the sequencer
            req = await self.seq_item_port.get_next_item()
            
            # Align to the falling edge before driving configuration signals
            await FallingEdge(self.dut.clk) 
            
            self.dut.rst_n.value = req.rst_n
            # Drive the configuration inputs FIRST to meet setup time requirements
            self.dut.FTW_start.value = req.FTW_start
            self.dut.FTW_step.value = req.FTW_step
            self.dut.cycles.value = req.cycles
            self.dut.rst_n.value = req.rst_n
            
            # Wait for the next falling edge to explicitly delay the enable signal
            await FallingEdge(self.dut.clk)
            self.dut.enable.value = req.enable 
          #  self.logger.info(f"Driving: {req.convert2string_stimulus()}")
            
            # Notify the sequencer that this item has been completely driven
            self.seq_item_port.item_done()