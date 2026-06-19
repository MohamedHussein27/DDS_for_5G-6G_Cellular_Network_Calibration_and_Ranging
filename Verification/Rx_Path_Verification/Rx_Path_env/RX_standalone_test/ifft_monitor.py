"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_monitor.py

    Description:
        The monitor is a passive component that observes the IFFT boundary.
        It reads the named internal wires of RX_TOP directly rather than the
        port handles of the u_rx_ifft submodule.

    WHY THE CHANGE:
        In Icarus Verilog (and most simulators via cocotb hierarchy), reading
        an INPUT port of a submodule handle (e.g. dut.u_rx_ifft.valid_in)
        always returns the port's own driven value — but cocotb cannot
        back-propagate the driving net (mult_valid) through the instance
        boundary. The result is that valid_in reads as 0 every cycle even
        when mult_valid is toggling, so the scoreboard sees 0 samples.

        The fix: give this monitor the top-level DUT handle (RX_TOP) and
        read the named internal wires directly:

            valid_in / in_real / in_imag  <- mult_valid / mult_re / mult_im
            valid_out / out_real / out_imag <- ifft_valid / ifft_re / ifft_im

        mult_re and mult_im are wires at RX_TOP scope; the >>>3 arithmetic
        shift in the RTL (.in_real(mult_re>>>3)) is not a named wire, so we
        replicate the >>3 here in Python after reading mult_re/mult_im.
"""

import cocotb
from cocotb.triggers import *
from pyuvm import *
from ifft_seq_item import ifft_item
from ifft_seq_item import *

def safe_int(val, default=0):
    try: return int(val)
    except ValueError: return default

def safe_signed(val, default=0):
    try: return val.signed_integer
    except ValueError: return default

#def arith_shift_right(val, n, bits=16):
#    """Replicate Verilog signed arithmetic right-shift (>>>n) in Python."""
    # sign-extend first
#    if val & (1 << (bits - 1)):
#        val -= (1 << bits)
#    return val >> n
def arith_shift_right(val, n):
    """Python's >> operator inherently performs arithmetic right-shift on signed ints."""
    return int(val) >> n



class ifft_monitor(uvm_monitor):
    def build_phase(self):
        # DUT handle must be the TOP-level RX_TOP handle (set in base_test.py
        # as ConfigDB "env.ifft_agt.*" -> self.dut, NOT self.dut.u_rx_ifft)
        self.dut = ConfigDB().get(self, "", "DUT")
        self.mon_port = uvm_analysis_port("mon_port", self)

    async def run_phase(self):
        await RisingEdge(self.dut.clk)

        while True:
            await RisingEdge(self.dut.clk)
            await ReadOnly()

            rsp_seq_item = ifft_item.create("rsp_seq_item")

            # ── Inputs to u_rx_ifft ─────────────────────────────────────
            # RTL: .valid_in(mult_valid)
            rsp_seq_item.valid_in = safe_int(self.dut.mult_valid.value)

            # RTL: .in_real(mult_re>>>3)  — read mult_re and shift in Python
            mult_re_raw = safe_signed(self.dut.mult_re.value)
            mult_im_raw = safe_signed(self.dut.mult_im.value)
            rsp_seq_item.in_real  = arith_shift_right(mult_re_raw, 3)
            rsp_seq_item.in_imag  = arith_shift_right(mult_im_raw, 3)

            # ── Outputs of u_rx_ifft ─────────────────────────────────────
            # RTL: .valid_out(ifft_valid), .out_real(ifft_re), .out_imag(ifft_im)
            rsp_seq_item.valid_out = safe_int(self.dut.ifft_valid.value)
            rsp_seq_item.out_real  = safe_signed(self.dut.ifft_re.value)
            rsp_seq_item.out_imag  = safe_signed(self.dut.ifft_im.value)
            # ── Control ──────────────────────────────────────────────────
            rsp_seq_item.rst_n = safe_int(self.dut.rst_n.value)

            self.logger.info(
                f"Monitor captured: rst_n={rsp_seq_item.rst_n}, "
                f"valid_in={rsp_seq_item.valid_in}, "
                f"in_real={rsp_seq_item.in_real}, in_imag={rsp_seq_item.in_imag}, "
                f"valid_out={rsp_seq_item.valid_out}, "
                f"out_real={rsp_seq_item.out_real}, out_imag={rsp_seq_item.out_imag}"
            )
            self.mon_port.write(rsp_seq_item)