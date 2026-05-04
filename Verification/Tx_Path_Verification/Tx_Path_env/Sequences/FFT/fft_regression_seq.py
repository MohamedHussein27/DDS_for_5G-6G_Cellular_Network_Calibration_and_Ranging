import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *
from pyuvm import uvm_sequence
from fft_seq_item import fft_item


from fft_dc_seq import SeqDc
from fft_impulse_seq import SeqImpulse
from fft_neg_dc_seq import SeqNegDc
from fft_fullscale_seq import SeqFullScale
from fft_ten_sample_seq import fft_ten_sample_seq
from fft_latency_seq import SeqLatency
from fft_tone_sweep_seq import SeqToneSweep
from fft_multiframe_seq import fft_two_frame_seq
from fft_coverage_closure_seq import FftCoverageClosureSeq
from max_neg_seq import FftUnscaledMinSeq
class SeqRegression(uvm_sequence):
    """
    Master Regression Sequence
    Description: Automatically runs every single FFT test back-to-back.
    """
    async def body(self):
        
        cocotb.log.info("=== STARTING MASTER REGRESSION SEQUENCE ===")

        # List of all your sequence classes and their display names
        sequence_list = [
            (SeqImpulse, "Impulse Test"),
            (SeqDc, "Positive DC Test"),
            (SeqNegDc, "Negative DC Test"),
            (SeqFullScale, "Nyquist Full Scale Test"),
            (fft_ten_sample_seq, "Ten Sample Zero-Pad Test"),
            (SeqLatency, "Random Latency/Throughput Test"),
            (fft_two_frame_seq, "Two Frame Continuous Test"),
            (SeqToneSweep, "5-Bin Tone Sweep Test"),
            (FftCoverageClosureSeq, "16-bit Corner Case Coverage Closure Test"),
            (FftUnscaledMinSeq, "Unscaled Minimum DC Test")
        ]

        for SeqClass, test_name in sequence_list:
            cocotb.log.info("==================================================")
            cocotb.log.info(f" EXECUTION: {test_name}")
            cocotb.log.info("==================================================")
            
            # Instantiate and start the sub-sequence
            seq = SeqClass("sub_seq")
            await seq.start(self.sequencer)

            cocotb.log.info(f"Flushing pipeline for {test_name}...")
            for _ in range(4500):
                await RisingEdge(cocotb.top.clk)

        cocotb.log.info("=== MASTER REGRESSION SEQUENCE COMPLETE ===")