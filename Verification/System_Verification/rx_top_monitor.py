import pyuvm
from pyuvm import *
import cocotb
from cocotb.triggers import RisingEdge, ReadOnly

# Assuming you named your sequence item class rx_item
from rx_top_seq_item import rx_top_item

class rx_monitor(uvm_monitor):
    def build_phase(self):
        self.mon_ap = uvm_analysis_port("mon_ap", self)
        self.dut_mon = ConfigDB().get(self, "", "DUT")
        
    async def run_phase(self):
        # Safely read 1-bit signals (handles 'X' and 'Z' at time 0)
        def read_1bit(sig):
            try:
                return int(sig.value)
            except ValueError:
                return 0
                
        # Safely read multi-bit signed signals
        def read_signed(sig):
            try:
                return sig.value.signed_integer
            except ValueError:
                return 0

        while True:
            await RisingEdge(self.dut_mon.clk)
            await ReadOnly()
            
            item = rx_top_item.create("item")
            
            # ==========================================
            # Capture Inputs (Stimulus)
            # ==========================================
            item.rst_n       = read_1bit(self.dut_mon.rst_n)
            
            # Channel Input Path (RESTORED TO MATCH NEW RTL)
            item.rx_valid_in = read_1bit(self.dut_mon.rx_valid_in)
            item.rx_in_re    = read_signed(self.dut_mon.rx_in_re)
            item.rx_in_im    = read_signed(self.dut_mon.rx_in_im)
            
            # Reference RAM Input Path
            item.ref_wr_en   = read_1bit(self.dut_mon.ref_wr_en)
            item.ref_wr_re   = read_signed(self.dut_mon.ref_wr_re)
            item.ref_wr_im   = read_signed(self.dut_mon.ref_wr_im)
            
            # ==========================================
            # Capture Outputs (Results)
            # ==========================================
            # OFDM Path
            item.ofdm_valid_out = read_1bit(self.dut_mon.ofdm_valid_out)
            item.ofdm_out_re    = read_signed(self.dut_mon.ofdm_out_re)
            item.ofdm_out_im    = read_signed(self.dut_mon.ofdm_out_im)
            
            # Radar Path
            item.radar_valid_out = read_1bit(self.dut_mon.radar_valid_out)
            item.radar_out_re    = read_signed(self.dut_mon.radar_out_re)
            item.radar_out_im    = read_signed(self.dut_mon.radar_out_im)
            
            self.mon_ap.write(item)