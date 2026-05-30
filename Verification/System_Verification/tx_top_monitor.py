import pyuvm
from pyuvm import *
from tx_top_seq_item import *
import cocotb
from cocotb.clock import Clock
from cocotb.triggers import *

class top_monitor(uvm_monitor):
    def build_phase(self):
        self.mon_ap = uvm_analysis_port("mon_ap", self)
        self.dut_mon = ConfigDB().get(self, "", "DUT")
        
    async def run_phase(self):
        await RisingEdge(self.dut_mon.clk)
        while True:
            await RisingEdge(self.dut_mon.clk)
            await ReadOnly()
            
            seq_item = tx_top_item.create("seq_item")
            
            # 1. Capture the Memory-Mapped Bus Inputs
            seq_item.rst_n   = int(self.dut_mon.rst_n.value)  
            seq_item.addr    = int(self.dut_mon.addr.value)
            seq_item.wr_en   = int(self.dut_mon.wr_en.value)
            seq_item.wr_data = int(self.dut_mon.wr_data.value)
            seq_item.rd_en   = int(self.dut_mon.rd_en.value)
            
            # 2. Capture the Bus Outputs (Read Data & Hardware Flags)
            # Using try/except in case rd_data or dds_ready_flag are 'X' or 'Z' during reset
            try:
                seq_item.rd_data = int(self.dut_mon.rd_data.value)
            except ValueError:
                seq_item.rd_data = 0
                
            try:
                seq_item.dds_ready_flag = int(self.dut_mon.dds_ready_flag.value)
            except ValueError:
                seq_item.dds_ready_flag = 0

            # 3. Capture the TX Datapath Outputs 
            if self.dut_mon.tx_valid.value == 1:
                seq_item.tx_valid    = 1
                seq_item.tx_out_real = int(self.dut_mon.tx_out_re.value)
                seq_item.tx_out_imag = int(self.dut_mon.tx_out_im.value)
            else:
                seq_item.tx_valid    = 0
                seq_item.tx_out_real = 0
                seq_item.tx_out_imag = 0
            
            # 4. Log and Send to Analysis Port
            self.logger.info(
                f"Monitoring BUS: rst={seq_item.rst_n}, addr={hex(seq_item.addr)}, we={seq_item.wr_en}, "
                f"wdata={hex(seq_item.wr_data)}, re={seq_item.rd_en} | "
                f"OUT: rdata={hex(seq_item.rd_data)}, ready={seq_item.dds_ready_flag}, "
                f"tx_valid={seq_item.tx_valid}"
            )
            
            self.mon_ap.write(seq_item)