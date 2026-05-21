"""
===============================================================================
Sponsor      : Analog Devices, Inc. (ADI)
Institution  : Faculty of Engineering, Ain Shams University
Project      : DDS for 5G/6G Cellular Network Calibration and Ranging

Module       : monitor.py

Description  :
    DDS monitor responsible for observing DUT activity and converting
    signal-level behavior into transaction-level sequence items.

    Main responsibilities:
        - Sample DUT signals every clock cycle
        - Capture DDS configuration and output data
        - Track chirp/sample indexing
        - Handle reset and invalid/X states safely
        - Broadcast collected transactions through the analysis port

===============================================================================
"""

import pyuvm
from pyuvm import *
from dds_seq_item import *

import cocotb
from cocotb.triggers import RisingEdge, ReadOnly


class dds_monitor(uvm_monitor):

    def build_phase(self):
        """
        Create monitor analysis port and retrieve DUT handle.
        """

        # Analysis port used to broadcast monitored transactions
        self.mon_ap = uvm_analysis_port("mon_ap", self)

        # Retrieve DUT handle from ConfigDB
        self.dut_mon = ConfigDB().get(self, "", "DUT")
        
    async def run_phase(self):
        """
        Main monitor sampling loop.
        """

        # Chirp/sample counter
        sample_index = 1

        while True:

            # Synchronize with DUT clock
            await RisingEdge(self.dut_mon.clk)

            # Wait for signal values to stabilize
            await ReadOnly()
            
            # -----------------------------------------------------------------
            # Create transaction object for current cycle
            # -----------------------------------------------------------------
            seq_item = dds_seq_item.create("seq_item")
            
            # -----------------------------------------------------------------
            # Capture reset signal
            # -----------------------------------------------------------------
            try:
                seq_item.rst_n = int(self.dut_mon.rst_n.value)

            except ValueError:
                # Handle X/Z states safely
                seq_item.rst_n = 0

            # -----------------------------------------------------------------
            # Capture control signals
            # -----------------------------------------------------------------
            try:
                seq_item.enable = int(self.dut_mon.enable.value)
                seq_item.valid_out = int(self.dut_mon.valid_out.value)

            except ValueError:
                seq_item.enable = 0
                seq_item.valid_out = 0

            # -----------------------------------------------------------------
            # Capture DDS output amplitude
            # -----------------------------------------------------------------
            try:
                val = int(self.dut_mon.final_amplitude.value)

                # Convert unsigned to signed 8-bit value
                if val > 127:
                    val -= 256

                seq_item.final_amplitude = val

            except ValueError:

                # Default handling during reset or invalid states
                seq_item.final_amplitude = (
                    0xFF if seq_item.rst_n == 0 else 0
                )

            # -----------------------------------------------------------------
            # Reset handling
            # -----------------------------------------------------------------
            if seq_item.rst_n == 0:

                # Reset chirp/sample tracking
                sample_index = 1
                seq_item.sample_index = 0

            else:

                # -------------------------------------------------------------
                # Capture DDS configuration parameters
                # -------------------------------------------------------------
                try:
                    seq_item.FTW_start = int(
                        self.dut_mon.FTW_start.value
                    )

                    seq_item.FTW_step = int(
                        self.dut_mon.FTW_step.value
                    )

                    seq_item.cycles = int(
                        self.dut_mon.cycles.value
                    )

                except ValueError:

                    # Handle X/Z values safely
                    seq_item.FTW_start = 0
                    seq_item.FTW_step = 0
                    seq_item.cycles = 0

                # Debug transaction print
                self.logger.debug(
                    f"Monitor captured: {seq_item.convert2string()}"
                )

                # Attach sample/chirp index
                seq_item.sample_index = sample_index

                # -------------------------------------------------------------
                # Update sample counter
                # -------------------------------------------------------------
                if seq_item.valid_out == 1:

                    # Increment while valid streaming continues
                    sample_index += 1

                else:

                    # Reset counter for next chirp/frame
                    sample_index = 1
                    
            # -----------------------------------------------------------------
            # Broadcast transaction to subscribers
            # -----------------------------------------------------------------
            self.mon_ap.write(seq_item)