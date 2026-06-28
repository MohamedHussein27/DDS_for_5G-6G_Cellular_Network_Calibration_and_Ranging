"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
 
    Module: tx_scoreboard.py
 
    Description:
        The TX scoreboard measures how closely the fixed-point RTL DUT
        matches the ideal floating-point TX signal using SQNR.
 
        FTW → Physical parameter recovery
        ──────────────────────────────────
            FTW_start = floor( f0 × 2^M / Fs )   →   f0 = FTW_start × Fs / 2^M
            FTW_step  = floor( B  × 2^M / (Fs×N) )  →   B  = FTW_step  × Fs × N / 2^M
 
        Reference model  (pure numpy built-ins, NO bit-true functions)
        ───────────────────────────────────────────────────────────────
        Stage 1 – Ideal DDS chirp
            n         = 1 … N  (1-indexed, same as RTL phase accumulator)
            phase[n]  = 2π × ( (f0/Fs)×n  +  (B/(Fs×N)) × n×(n−1)/2 )
            s[n]      = A × sin( phase[n] )     A = 127/256
 
        Stage 2 – np.fft.fft(s)
 
        Stage 3 – Ideal MUX
            OFDM integers from DUT ROM are in Q8.8 → de-convert to float
            first (divide by 256), then place in the MUX frame as floats.
 
            Bin layout (0-based, verified against standalone reference run):
              Bins 0    – 2047 : OFDM symbols  (de-converted to float)
              Bins 2048 – 2429 : zeros          (382 bins, index 2429 is zero)
              Bins 2430 – 4095 : FFT / 128.0    (1666 bins)
 
            NOTE: the +1 shift (2430 instead of 2429) was confirmed by
            comparing _build_reference() output against the standalone
            reference script: index 2429 is zero and the first non-zero
            radar bin appears at index 2430.
 
        Stage 4 – np.fft.ifft(mux_frame) × N
            ×N because the RTL IFFT does NOT divide by N.
 
        Stage 5 – Convert IFFT output to Q8.8
            ref_int = clip( round( x × 256 ), −32768, 32767 )
 
        Comparison metric (SQNR)
        ────────────────────────
            SQNR = 10 × log10( Σ|ref_int|² / Σ|dut_int − ref_int|² )
 
        Pass/fail thresholds:
            SQNR_PASS_DB = 30 dB
            SQNR_WARN_DB = 25 dB
"""
 
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import logging
import numpy as np
 
from top_seq_item import *
 
# ─────────────────────────────────────────────────────────────────────────────
# System constants
# ─────────────────────────────────────────────────────────────────────────────
Fs           = 491_520_000.0    # sampling frequency (Hz)
M            = 32               # phase-accumulator width (bits)
N            = 4096             # frame size (samples)
A_DDS        = 127.0 / 256.0   # DDS amplitude
theta        = 0.0              # initial phase (rad)
 
Q8_8_SCALE   = 256.0           # 2^FL_Q8_8 = 2^8
 
# MUX layout — verified against standalone reference run:
#   index 2429 is zero, first radar bin is at 2430
N_OFDM       = 2048
RADAR_START  = 2429             # 2429 is zero; radar begins at 2430
N_RADAR      = N - RADAR_START  # 4096 - 2430 = 1666
RADAR_SHIFT  = 128.0            # floating-point equivalent of RTL >>7
 
# SQNR thresholds
SQNR_PASS_DB = 30.0
SQNR_WARN_DB = 25.0
 
 
# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def _sign16(v):
    """Sign-extend a raw 16-bit integer from the DUT."""
    v = int(v) & 0xFFFF
    return (v - 0x10000) if (v & 0x8000) else v
 
 
def _q8_8_to_float(int_array):
    """
    De-convert Q8.8 signed integers back to floating-point.
    Mirrors convert_hex_q8_8_to_float() exactly:
        float = sign16(raw_int) / 256.0
    The _sign16 step handles two's complement for values ≥ 32768.
    """
    arr = np.array(int_array, dtype=np.int64)
    # Apply two's complement (same as the deconverter: subtract 65536 if ≥ 32768)
    arr[arr >= 32768] -= 65536
    return arr / Q8_8_SCALE
 
 
def _to_q8_8(arr_float):
    """
    Convert floating-point array to Q8.8 signed integers.
    Mirrors convert_complex_float_to_q8_8() exactly:
        q = clip( round( x × 256 ), −32768, 32767 )
    """
    return np.clip(
        np.round(arr_float * Q8_8_SCALE), -32768, 32767
    ).astype(np.int64)
 
 
def _recover_f0_B(FTW_start, FTW_step):
    """
    Invert the hand-written FTW equations:
        f0 = FTW_start × Fs / 2^M
        B  = FTW_step  × Fs × N / 2^M
    """
    scale = float(1 << M)      # 2^32
    f0    = FTW_start * Fs / scale
    B     = FTW_step  * Fs * N / scale
    return f0, B
 
 
# ─────────────────────────────────────────────────────────────────────────────
# Ideal floating-point reference pipeline  (numpy built-ins only)
# ─────────────────────────────────────────────────────────────────────────────
def _build_reference(FTW_start, FTW_step, cycles,  ofdm_re_ints, ofdm_im_ints):
    """
    Recover f0/B from FTW values, compute the ideal TX output,
    export ALL intermediate stages to text files for debugging,
    and return Q8.8 integer reference vectors.
    """

    # ── Recover physical parameters ───────────────────────────────────────────
    f0, B = _recover_f0_B(FTW_start, FTW_step)
    print(f"[REF] FTW_start={FTW_start}  FTW_step={FTW_step}")
    print(f"[REF] Recovered f0={f0/1e6:.6f} MHz   B={B/1e6:.6f} MHz")

    # ── Stage 1: Ideal DDS chirp ──────────────────────────────────────────────
    n         = np.arange(1, cycles + 1, dtype=np.float64)  
    phase_rad = 2.0 * np.pi * (
        (f0 / Fs) * n
        + (B / (Fs * N)) * (n * (n - 1) / 2.0)
    )
    s = A_DDS * np.sin(phase_rad + theta)
    
    # EXPORT DDS
    np.savetxt("py_dbg_stage1_dds_out.txt", s, fmt='%.6f')

    # ── Stage 2: FFT ──────────────────────────────────────────────────────────
    fft_out = np.fft.fft(s, n=N)

    # EXPORT FFT (Split into Real and Imaginary floats)
    np.savetxt("py_dbg_stage2_fft_re.txt", np.real(fft_out), fmt='%.6f')
    np.savetxt("py_dbg_stage2_fft_im.txt", np.imag(fft_out), fmt='%.6f')

    # ── Stage 3: Ideal MUX ────────────────────────────────────────────────────
    ofdm_float_re = _q8_8_to_float(ofdm_re_ints)
    ofdm_float_im = _q8_8_to_float(ofdm_im_ints)

    # Shift values: fill index 1 to the end with values from index 0 to second-to-last
    ofdm_float_re[1:] = ofdm_float_re[:-1]
    ofdm_float_im[1:] = ofdm_float_im[:-1]

    ofdm_complex  = ofdm_float_re + 1j * ofdm_float_im

    mux_frame = np.zeros(N, dtype=np.complex128)

    # State 1: OFDM
    mux_frame[:N_OFDM] = ofdm_complex[:N_OFDM]

    # Enforce Guardband 1
    bin_res = Fs / N
    gb1_start = int(200e3 / bin_res)
    gb1_end   = int(210e3 / bin_res)
    mux_frame[gb1_start:gb1_end+1] = 0

    # State 3: Radar/Chirp
    mux_frame[RADAR_START:] = fft_out[:N_RADAR] / RADAR_SHIFT

    # EXPORT MUX (Split into Real and Imaginary floats)
    np.savetxt("py_dbg_stage3_mux_re.txt", np.real(mux_frame), fmt='%.6f')
    np.savetxt("py_dbg_stage3_mux_im.txt", np.imag(mux_frame), fmt='%.6f')

    # ── Stage 4: IFFT ─────────────────────────────────────────────────────────
    ifft_unscaled = np.fft.ifft(mux_frame) * N

    ifft_re = np.real(ifft_unscaled)
    ifft_im = np.imag(ifft_unscaled)

    # EXPORT IFFT (Unscaled floats)
    np.savetxt("py_dbg_stage4_ifft_re.txt", ifft_re, fmt='%.6f')
    np.savetxt("py_dbg_stage4_ifft_im.txt", ifft_im, fmt='%.6f')

    # ── Stage 5: Convert IFFT output to Q8.8 ─────────────────────────────────
    ref_int_re = _to_q8_8(ifft_re)
    ref_int_im = _to_q8_8(ifft_im)

    # ── Stage 6: Export Final Q8.8 Reference Vectors ─────────────────────────
    np.savetxt("py_ref_ifft_out_q8_8_re.txt", ref_int_re, fmt='%d')
    np.savetxt("py_ref_ifft_out_q8_8_im.txt", ref_int_im, fmt='%d')

    print("[REF] All intermediate stages successfully dumped to .txt files.")

    return ref_int_re, ref_int_im


# ─────────────────────────────────────────────────────────────────────────────
# Scoreboard
# ─────────────────────────────────────────────────────────────────────────────
class top_scoreboard_sqnr(uvm_scoreboard):
    """
    Scoreboard for the full TX path (TX_TOP).

    Inputs captured by the monitor in top_item:
        rst_n, dds_ready_flag

    Inputs captured by the monitor in dds_item:
        enable, FTW_start, FTW_step, cycles

    Outputs captured by the monitor in top_item:
        tx_valid, tx_out_real, tx_out_imag
    """

    def __init__(self, name, parent):
        super().__init__(name, parent)

    # ── build_phase ───────────────────────────────────────────────────────────
    def build_phase(self):
        # Analysis export/FIFO for top_item transactions
        self.sb_export  = uvm_analysis_export("sb_export", self)
        self.sb_fifo    = uvm_tlm_analysis_fifo("sb_fifo", self)
        self.sb_export  = self.sb_fifo.analysis_export

        # Analysis export/FIFO for dds_item transactions
        self.dds_export = uvm_analysis_export("dds_export", self)
        self.dds_fifo   = uvm_tlm_analysis_fifo("dds_fifo", self)
        self.dds_export = self.dds_fifo.analysis_export

        # DUT handle for OFDM ROM access
        self.dut = ConfigDB().get(self, "", "DUT")

        # Counters
        self.reset_cycles  = 0

        # Reference vectors (computed once on first dds_enable)
        self._ref_int_re   = None
        self._ref_int_im   = None
        self._golden_ready = False
        self._out_idx      = 0

        # SQNR accumulators
        self._sig_pow_acc   = 0.0   # Σ ( ref_int_re² + ref_int_im² )
        self._noise_pow_acc = 0.0   # Σ ( (dut_re−ref_re)² + (dut_im−ref_im)² )

    # ── run_phase ─────────────────────────────────────────────────────────────
    async def run_phase(self):
        while True:
            item     = await self.sb_fifo.get()
            dds_item = await self.dds_fifo.get()
            self._process(item, dds_item)

    # ── _process ──────────────────────────────────────────────────────────────
    def _process(self, item, dds_item):

        global i

        # ── Reset ─────────────────────────────────────────────────────────────
        if not item.rst_n:
            self.reset_cycles      += 1
            self._golden_ready      = False
            self._ref_int_re        = None
            self._ref_int_im        = None
            self._out_idx           = 0
            self._sig_pow_acc       = 0.0
            self._noise_pow_acc     = 0.0
            self.logger.info("Reset detected – scoreboard state flushed.")
            i = 0
            return

        # ── Build ideal reference once on first dds_enable ────────────────────
        if dds_item.enable and not self._golden_ready:
            self.logger.info(
                f"Building ideal floating-point reference: "
                f"FTW_start={dds_item.FTW_start}, "
                f"FTW_step={dds_item.FTW_step}, "
                f"cycles={dds_item.cycles}"
            )

            # Read OFDM ROM directly from the DUT
            ofdm_re_list = []
            ofdm_im_list = []
            rom_re    = self.dut.u_ofdm_rom.rom_real
            rom_im    = self.dut.u_ofdm_rom.rom_imag
            rom_depth = len(rom_re)

            print(f"Reading OFDM ROM (depth={rom_depth})...")
            for k in range(rom_depth):
                try:
                    val_re = _sign16(int(rom_re[k].value))
                    val_im = _sign16(int(rom_im[k].value))
                except ValueError:
                    val_re, val_im = 0, 0
                ofdm_re_list.append(val_re)
                ofdm_im_list.append(val_im)

            # Run the ideal pipeline
            self._ref_int_re, self._ref_int_im = _build_reference(
                dds_item.FTW_start,
                dds_item.FTW_step,
                dds_item.cycles,
                ofdm_re_list,
                ofdm_im_list,
            )

            self._golden_ready = True
            self._out_idx      = 0

        # ── Heartbeat / debug prints ───────────────────────────────────────────
        if (i > 4090 and i < 4100) or (i > 12275 and i < 12285):
            self.logger.info(
                f"[i={i}] DUT: re={_sign16(item.tx_out_real)}, "
                f"im={_sign16(item.tx_out_imag)}, "
                f"dds_ready={item.dds_ready_flag}, tx_valid={item.tx_valid}"
            )

        if i % 1000 == 0:
            self.logger.info(f"--- SIMULATION HEARTBEAT: Processing cycle {i} ---")

        # ── Collect DUT output and accumulate SQNR ────────────────────────────
        if item.tx_valid:
            if not self._golden_ready:
                self.logger.warning(
                    "tx_valid asserted before reference is ready – skipping."
                )
            elif self._out_idx >= N:
                self.logger.warning(
                    f"tx_valid beyond frame size N={N} "
                    f"(out_idx={self._out_idx}) – skipping."
                )
            else:
                dut_re = _sign16(item.tx_out_real)
                dut_im = _sign16(item.tx_out_imag)
                ref_re = int(self._ref_int_re[self._out_idx])
                ref_im = int(self._ref_int_im[self._out_idx])

                err_re = dut_re - ref_re
                err_im = dut_im - ref_im

                # Debug at frame boundaries
                if i > 4090 and i < 4100:
                    self.logger.info(
                        f"Sample #{self._out_idx}: "
                        f"DUT=({dut_re}, {dut_im})  "
                        f"REF=({ref_re}, {ref_im})  "
                        f"ERR=({err_re}, {err_im})"
                    )

                # Accumulate signal and noise power
                self._sig_pow_acc   += ref_re**2  + ref_im**2
                self._noise_pow_acc += err_re**2  + err_im**2

                self._out_idx += 1

        i += 1

    # ── report_phase ──────────────────────────────────────────────────────────
    def report_phase(self):

        # ── Compute SQNR ──────────────────────────────────────────────────────
        if self._noise_pow_acc > 0:
            sqnr_db = 10.0 * np.log10(self._sig_pow_acc / self._noise_pow_acc)
        else:
            sqnr_db = float("inf")

        sqnr_str = f"{sqnr_db:.2f} dB" if sqnr_db != float("inf") else "∞  (perfect match)"

        # ── Pass / Warn / Fail decision ───────────────────────────────────────
        if sqnr_db >= SQNR_PASS_DB:
            verdict = "PASS ✓"
        elif sqnr_db >= SQNR_WARN_DB:
            verdict = "WARN ⚠"
        else:
            verdict = "FAIL ✗"

        # ── Report ────────────────────────────────────────────────────────────
        self.logger.info("╔══════════════════════════════════════════════════════╗")
        self.logger.info("║              TX SCOREBOARD REPORT                    ║")
        self.logger.info("╠══════════════════════════════════════════════════════╣")
        self.logger.info(f"║  Reset cycles detected        : {self.reset_cycles:<22} ║")
        self.logger.info(f"║  Samples compared             : {self._out_idx:<22} ║")
        self.logger.info("║ ───────────────────────────────────────────────────  ║")
        self.logger.info(f"║  Signal power  Σ|ref|²        : {self._sig_pow_acc:<22.2f} ║")
        self.logger.info(f"║  Noise power   Σ|err|²        : {self._noise_pow_acc:<22.2f} ║")
        self.logger.info(f"║  SQNR                         : {sqnr_str:<22} ║")
        self.logger.info(f"║  SQNR threshold (PASS)        : {SQNR_PASS_DB:<22} ║")
        self.logger.info(f"║  SQNR threshold (WARN)        : {SQNR_WARN_DB:<22} ║")
        self.logger.info("║ ───────────────────────────────────────────────────  ║")
        self.logger.info(f"║  OVERALL RESULT               : {verdict:<22} ║")
        self.logger.info("╚══════════════════════════════════════════════════════╝")

        # Hard-fail or warn
        if sqnr_db < SQNR_WARN_DB:
            self.logger.critical(
                f"SQNR = {sqnr_db:.2f} dB is below the FAIL threshold "
                f"of {SQNR_WARN_DB} dB."
            )
        elif sqnr_db < SQNR_PASS_DB:
            self.logger.warning(
                f"SQNR = {sqnr_db:.2f} dB is in the warning zone "
                f"({SQNR_WARN_DB}–{SQNR_PASS_DB} dB)."
            )