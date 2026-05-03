"""
    Module: rx_monitor.py

    Description:
        The RX monitor passively observes all activity on the RX_TOP boundaries.
        It samples inputs and outputs synchronously on the rising edge of the clock
        (after ReadOnly to avoid race conditions) and broadcasts the complete 
        transaction via an analysis port to the Scoreboard and Subscribers.
"""

import pyuvm
from pyuvm import *
import cocotb
from cocotb.triggers import RisingEdge, ReadOnly

# Assuming you named your sequence item class rx_item
from top_seq_item import rx_item

class rx_monitor(uvm_monitor):
    def build_phase(self):
        # Create the analysis port to broadcast transactions
        self.mon_ap = uvm_analysis_port("mon_ap", self)
        # Fetch the DUT from the Configuration Database
        self.dut_mon = ConfigDB().get(self, "", "DUT")
        
    async def run_phase(self):
        while True:
            # 1. Synchronize to the active clock edge
            await RisingEdge(self.dut_mon.clk)
            
            # 2. Wait until all signals have settled in this delta cycle
            await ReadOnly()
            
            # 3. Create a new transaction item
            # Note: Do not use seq_item = seq_item.create(), use the class name
            item = rx_item("mon_item")
            
            # ==========================================
            # 4. Capture Inputs (Stimulus)
            # ==========================================
            # Global
            item.rst_n = self.dut_mon.rst_n.value.integer
            
            # Channel Input Path
            item.rx_valid_in = self.dut_mon.rx_valid_in.value.integer
            item.rx_in_re = self.dut_mon.rx_in_re.value.signed_integer
            item.rx_in_im = self.dut_mon.rx_in_im.value.signed_integer
            
            # Reference RAM Input Path
            item.ref_wr_en = self.dut_mon.ref_wr_en.value.integer
            item.ref_wr_re = self.dut_mon.ref_wr_re.value.signed_integer
            item.ref_wr_im = self.dut_mon.ref_wr_im.value.signed_integer
            
            # ==========================================
            # 5. Capture Outputs (Results)
            # ==========================================
            # OFDM Path
            item.ofdm_valid_out = self.dut_mon.ofdm_valid_out.value.integer
            item.ofdm_out_re = self.dut_mon.ofdm_out_re.value.signed_integer
            item.ofdm_out_im = self.dut_mon.ofdm_out_im.value.signed_integer
            
            # Radar Path
            item.radar_valid_out = self.dut_mon.radar_valid_out.value.integer
            item.radar_out_re = self.dut_mon.radar_out_re.value.signed_integer
            item.radar_out_im = self.dut_mon.radar_out_im.value.signed_integer
            
            # ==========================================
            # 6. Capture Internal Probes (Optional/If exposed)
            # ==========================================
            # If these are exposed to the top level for white-box monitoring:
            # item.fft_valid = self.dut_mon.fft_valid.value.integer
            # item.radar_split_valid = self.dut_mon.radar_split_valid.value.integer
            # item.mult_valid = self.dut_mon.mult_valid.value.integer
            # item.ifft_valid = self.dut_mon.ifft_valid.value.integer

            # 7. Broadcast the completed transaction to connected components
            self.mon_ap.write(item)