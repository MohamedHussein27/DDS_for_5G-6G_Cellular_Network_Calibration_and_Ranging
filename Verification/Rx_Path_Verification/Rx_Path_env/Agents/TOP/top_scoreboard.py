"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: tx_scoreboard.py

    Description:
        The TX-path scoreboard receives paired transactions from the monitor.
        The full TX pipeline is:

            DDS (chirp time-domain)
              → FFT (chirp spectrum, Q8.8 after reinterpret-cast)
              → Bit-Reversal (bit-reversed FFT frame)
              → ISAC MUX  (2048 OFDM | 381 zeros | 1667 radar bins)
              → IFFT (final TX time-domain output, Q8.8)

        Because every stage has its own latency, the monitor captures all
        intermediate signals on a single clock edge and packs them into one
        TxItem.  The scoreboard checks each stage independently:

          Stage 1 – DDS amplitude: range sanity check (signed 8-bit)
          Stage 2 – FFT output   : compared against numpy FFT + Q8.8 cast
          Stage 3 – Bit-reversal : output[k] == fft_frame[bit_rev(k)]
          Stage 4 – MUX frame    : structural check (OFDM / zeros / radar)
          Stage 5 – TX output    : compared against golden IFFT, SQNR ≥ 30 dB

        Word-length constants mirror Tx_path_fixed.py:
            WL     = 16    (all pipeline stages)
            FL_FFT =  8    (Q8.8, after reinterpret-cast with scale=64→log2=6)
            FL_TX  =  8    (TX IFFT output, same WL/FL as FFT input)
            TOLERANCE = 4  (±LSBs, loosest needed by the IFFT output stage)
"""

import numpy as np
from collections import deque
from pyuvm import *

from tx_item import TxItem          # your monitor transaction class

# ─────────────────────────────────────────────────────────────────────────────
# RTL / golden-model parameters  (mirror Tx_path_fixed.py exactly)
# ─────────────────────────────────────────────────────────────────────────────
WL          = 16    # word length for all pipeline stages
FL_FFT      =  8    # fractional bits after FFT reinterpret-cast  (Q8.8)
FL_TX       =  8    # fractional bits at the TX final output       (Q8.8)
N           = 4096  # FFT / IFFT / MUX frame size
ADDR_W      =  12   # log2(N) – used for bit-reversal address

# MUX output frame segment boundaries (0-based sample index)
SEG_OFDM_START  =    0
SEG_OFDM_END    = 2047   # 2048 OFDM bins
SEG_ZERO_START  = 2048
SEG_ZERO_END    = 2428   # 381 guard/padding zeros
SEG_RADAR_START = 2429
SEG_RADAR_END   = 4095   # 1667 radar chirp bins

TOLERANCE_DDS    = 0     # DDS is a range check only – no numeric tolerance
TOLERANCE_FFT    = 2     # ±2 LSBs in Q8.8
TOLERANCE_BITREV = 0     # bit-reversal is lossless – exact match
TOLERANCE_MUX    = 2     # ±2 LSBs in Q8.8
TOLERANCE_TX     = 4     # ±4 LSBs in Q8.8  (accumulated IFFT quantisation)

SQNR_THRESHOLD_DB = 30.0  # minimum acceptable TX output SQNR


# ─────────────────────────────────────────────────────────────────────────────
# Fixed-point helpers
# ─────────────────────────────────────────────────────────────────────────────
def _sign16(v: int) -> int:
    """Sign-extend a raw 16-bit integer value coming from the DUT."""
    v = int(v) & 0xFFFF
    return v - 0x10000 if (v & 0x8000) else v


def _sign8(v: int) -> int:
    """Sign-extend a raw 8-bit integer value (DDS amplitude)."""
    v = int(v) & 0xFF
    return v - 0x100 if (v & 0x80) else v


def _raw_to_float(v: int, frac_bits: int) -> float:
    """Convert a raw signed fixed-point integer to a Python float."""
    return _sign16(v) / (1 << frac_bits)


def _float_to_raw(f: float, frac_bits: int) -> int:
    """
    Convert a float to a 16-bit signed fixed-point integer.
    Uses convergent (banker's) rounding and saturates to [-32768, 32767].
    """
    raw = int(np.round(f * (1 << frac_bits)))
    return max(-(1 << (WL - 1)), min((1 << (WL - 1)) - 1, raw))


def _bit_reverse(idx: int, width: int = ADDR_W) -> int:
    """Reverse the 'width' LSBs of idx."""
    result = 0
    for _ in range(width):
        result = (result << 1) | (idx & 1)
        idx >>= 1
    return result


def _quantize_q8_8(x: np.ndarray) -> np.ndarray:
    """
    Quantize a complex float array to Q8.8 fixed-point and back to float.
    Convergent rounding + saturation – mirrors fi_cast in Tx_path_fixed.py.
    """
    scale   = 1 << FL_FFT         # 256
    max_val = ((1 << (WL - 1)) - 1) / scale
    min_val = -(1 << (WL - 1))    / scale

    def _q(v):
        rounded = np.round(v * scale) / scale
        return np.clip(rounded, min_val, max_val)

    return _q(x.real) + 1j * _q(x.imag)


# ─────────────────────────────────────────────────────────────────────────────
# Per-stage golden models
# ─────────────────────────────────────────────────────────────────────────────

class _DdsGolden:
    """
    Stage 1: DDS amplitude is a signed 8-bit sine sample.
    We perform a range check only (no LUT rebuild here).
    """
    @staticmethod
    def check(raw_amplitude: int) -> bool:
        """Return True if the 8-bit signed amplitude is in [-128, 127]."""
        amp = _sign8(raw_amplitude)
        return -128 <= amp <= 127


class _FftGolden:
    """
    Stage 2: Accumulate N DDS samples, run numpy FFT, apply the
    reinterpret-cast from Tx_path_fixed.py, and yield Q8.8 bins.

    Reinterpret-cast logic (mirrors the Python golden model):
        log2_scale = 6  (scale = 64)
        dst_FL     = FL_FFT + log2_scale = 8 + 6 = 14   ← src FL after FFT
        The raw integer is unchanged; only the FL label changes.
        Numerically: output_float = raw_integer / 2^14 → re-expressed as Q8.8
        i.e. value = (raw_int / 2^14) which we then store back as Q8.8.
    """
    LOG2_SCALE = 6       # log2(64)
    FL_SRC     = FL_FFT  # FL label before cast  (8 after fft, before cast)
    FL_DST     = FL_FFT + LOG2_SCALE   # = 14  (label after cast)

    def __init__(self):
        self._dds_samples = []   # raw signed 8-bit integers

    def push_dds(self, raw_amplitude: int):
        """Accept one DDS sample (raw signed 8-bit)."""
        self._dds_samples.append(_sign8(raw_amplitude))

    def ready(self) -> bool:
        return len(self._dds_samples) == N

    def compute(self) -> list:
        """
        Compute the chirp FFT + reinterpret-cast.
        Returns a list of N (re_raw_int, im_raw_int) Q8.8 pairs.
        """
        # Sign-extended 8-bit → Q1.0 float → sign-extended to 16-bit Q8.8
        # RTL does:  fft_in_re = { {8{dds_amplitude[7]}}, dds_amplitude }
        # i.e. sign-extend 8-bit → 16-bit.  FL stays at 0 (integer signal).
        dds_float = np.array(self._dds_samples, dtype=float)

        # FFT (numpy, equivalent to hardware radix-2 DIF)
        spectrum = np.fft.fft(dds_float)

        # Scale by 1/N to match the hardware's non-normalized convention
        # then re-express at FL_DST and round back to FL_FFT.
        # In practice the hardware just shifts bits; we replicate the
        # reinterpretcast numerically: divide by 2^(FL_DST - FL_SRC).
        cast_scale = 2 ** (self.FL_DST - self.FL_SRC)   # 2^6 = 64
        spectrum_cast = spectrum / cast_scale

        # Quantize to Q8.8
        spectrum_q = _quantize_q8_8(spectrum_cast)

        self._dds_samples.clear()

        result = []
        for v in spectrum_q:
            result.append((_float_to_raw(v.real, FL_FFT),
                           _float_to_raw(v.imag, FL_FFT)))
        return result


class _BitRevGolden:
    """
    Stage 3: Collect one full FFT frame, output in bit-reversed order.
    bit_rev_out[k] must equal fft_out[bit_rev(k)].
    """
    def __init__(self):
        self._fft_frame = []   # list of (re_raw, im_raw) – filled from FFT checker

    def load_fft_frame(self, fft_pairs: list):
        """Load a full N-sample FFT frame for later bit-reversal checking."""
        self._fft_frame = list(fft_pairs)

    def get_expected(self, output_idx: int):
        """Return the expected (re_raw, im_raw) for output sample output_idx."""
        if not self._fft_frame:
            return None
        src_idx = _bit_reverse(output_idx, ADDR_W)
        return self._fft_frame[src_idx]


class _MuxGolden:
    """
    Stage 4: Knows the three segments of a MUX output frame.
    Checks OFDM data, zero padding, and radar bins separately.
    """
    def __init__(self):
        self._ofdm_frame  = []    # 2048 Q8.8 (re_raw, im_raw) pairs
        self._radar_frame = []    # 1667 Q8.8 (re_raw, im_raw) pairs (after >>>7)

    def load_ofdm(self, pairs: list):
        """Load 2048 OFDM golden pairs (Q8.8 raw integers)."""
        self._ofdm_frame = list(pairs)

    def load_radar(self, pairs: list):
        """Load 1667 radar golden pairs (bit_rev output already >>>7, Q8.8)."""
        self._radar_frame = list(pairs)

    def get_expected(self, sample_idx: int):
        """
        Return (re_raw, im_raw, segment_name) for a given MUX sample index.
        Returns (0, 0, 'ZERO') for the guard band.
        Returns None if idx is out of range.
        """
        if SEG_OFDM_START <= sample_idx <= SEG_OFDM_END:
            if self._ofdm_frame:
                return (*self._ofdm_frame[sample_idx], 'OFDM')
            return None
        elif SEG_ZERO_START <= sample_idx <= SEG_ZERO_END:
            return (0, 0, 'ZERO')
        elif SEG_RADAR_START <= sample_idx <= SEG_RADAR_END:
            ridx = sample_idx - SEG_RADAR_START
            if self._radar_frame:
                return (*self._radar_frame[ridx], 'RADAR')
            return None
        return None


class _TxGolden:
    """
    Stage 5: Accumulate N MUX output samples, run numpy IFFT,
    quantize to Q8.8, and yield the expected TX time-domain samples.
    """
    def __init__(self):
        self._mux_samples = []   # complex floats (from Q8.8 raw pairs)
        self._out = deque()

    def push_mux(self, re_raw: int, im_raw: int):
        """Accept one MUX output sample."""
        re = _sign16(re_raw) / (1 << FL_FFT)
        im = _sign16(im_raw) / (1 << FL_FFT)
        self._mux_samples.append(complex(re, im))

    def ready(self) -> bool:
        return len(self._mux_samples) == N

    def compute(self):
        """Run the golden IFFT and buffer all N output samples."""
        x_in = np.array(self._mux_samples)

        # Hardware IFFT = numpy ifft * N  (hardware does not divide by N)
        ifft_out = np.fft.ifft(x_in) * N

        # Quantize to Q8.8
        ifft_q = _quantize_q8_8(ifft_out)

        self._mux_samples.clear()
        for v in ifft_q:
            self._out.append((_float_to_raw(v.real, FL_TX),
                              _float_to_raw(v.imag, FL_TX)))

    def pop(self):
        """Return the next (re_raw, im_raw) pair or None if not ready."""
        return self._out.popleft() if self._out else None


# ─────────────────────────────────────────────────────────────────────────────
# TX Scoreboard
# ─────────────────────────────────────────────────────────────────────────────
class TxScoreboard(uvm_scoreboard):
    """
    Scoreboard for the full TX path (TX_TOP).

    Expected TxItem fields (set by your monitor):
        rst_n           – reset level
        # DDS
        dds_valid       – DDS output valid
        dds_amplitude   – raw 8-bit signed amplitude
        # FFT
        fft_valid       – FFT output valid
        fft_re, fft_im  – raw 16-bit signed Q8.8
        # Bit-reversal (exposed debug ports)
        bit_rev_valid   – bit-reversal valid
        bit_rev_re      – raw 16-bit signed Q8.8
        bit_rev_im      – raw 16-bit signed Q8.8
        # MUX
        mux_valid       – MUX output valid
        mux_re, mux_im  – raw 16-bit signed Q8.8
        # TX final output
        tx_valid        – TX IFFT output valid
        tx_out_re       – raw 16-bit signed Q8.8
        tx_out_im       – raw 16-bit signed Q8.8
    """

    def __init__(self, name, parent):
        super().__init__(name, parent)

    # ── build_phase ───────────────────────────────────────────────────────────
    def build_phase(self):
        self.sb_export = uvm_analysis_export("sb_export", self)
        self.sb_fifo   = uvm_tlm_analysis_fifo("sb_fifo", self)
        self.sb_export = self.sb_fifo.analysis_export

        # ── Counters (one correct/error pair per stage) ───────────────────────
        self.correct = {"DDS": 0, "FFT_re": 0, "FFT_im": 0,
                        "BITREV_re": 0, "BITREV_im": 0,
                        "MUX_re": 0, "MUX_im": 0,
                        "TX_re": 0, "TX_im": 0}
        self.errors  = {"DDS": 0, "FFT_re": 0, "FFT_im": 0,
                        "BITREV_re": 0, "BITREV_im": 0,
                        "MUX_re": 0, "MUX_im": 0,
                        "TX_re": 0, "TX_im": 0}
        self.reset_cycles = 0

        # ── Golden model instances ────────────────────────────────────────────
        self._dds_gold    = _DdsGolden()
        self._fft_gold    = _FftGolden()
        self._bitrev_gold = _BitRevGolden()
        self._mux_gold    = _MuxGolden()
        self._tx_gold     = _TxGolden()

        # Internal state
        self._fft_golden_frame = []   # last computed FFT golden pairs
        self._fft_out_idx      = 0    # which FFT output bin we are at
        self._bitrev_out_idx   = 0    # which bit-rev output sample we are at
        self._mux_out_idx      = 0    # position within current MUX frame
        self._bitrev_frame     = []   # bit-rev output collected for MUX radar ref

        # SQNR accumulators (TX output stage)
        self._sqnr_sig_pow   = 0.0
        self._sqnr_noise_pow = 0.0

        # Sample counter (for debug prints – mirrors the IFFT scoreboard style)
        self._i = 0

    # ── run_phase ─────────────────────────────────────────────────────────────
    async def run_phase(self):
        while True:
            item = await self.sb_fifo.get()
            self._process(item)

    # ── _process ──────────────────────────────────────────────────────────────
    def _process(self, item: TxItem):

        # ── Reset: flush all golden model state ──────────────────────────────
        if not item.rst_n:
            self.reset_cycles += 1
            self._fft_gold        = _FftGolden()
            self._bitrev_gold     = _BitRevGolden()
            self._mux_gold        = _MuxGolden()
            self._tx_gold         = _TxGolden()
            self._fft_golden_frame.clear()
            self._fft_out_idx     = 0
            self._bitrev_out_idx  = 0
            self._mux_out_idx     = 0
            self._bitrev_frame.clear()
            self._i               = 0
            self.logger.info("Reset detected – all golden models flushed.")
            return

        # ── Debug print (mirrors the IFFT scoreboard style) ──────────────────
        if self._i == 0:
            self.logger.info(
                f"[i=0] dds_valid={item.dds_valid} fft_valid={item.fft_valid} "
                f"mux_valid={item.mux_valid} tx_valid={item.tx_valid}"
            )
        if 4090 < self._i < 4100:
            self.logger.info(
                f"[i={self._i}] TX out: re={_sign16(item.tx_out_re)} "
                f"im={_sign16(item.tx_out_im)} tx_valid={item.tx_valid}"
            )
        self._i += 1

        # ═══════════════════════════════════════════════════════════════════
        # STAGE 1 – DDS amplitude range check
        # ═══════════════════════════════════════════════════════════════════
        if item.dds_valid:
            ok = self._dds_gold.check(item.dds_amplitude)
            if ok:
                self.correct["DDS"] += 1
            else:
                self.errors["DDS"] += 1
                self.logger.error(
                    f"MISMATCH DDS amplitude "
                    f"(error #{self.errors['DDS']}) | "
                    f"DUT={_sign8(item.dds_amplitude)} out of signed 8-bit range"
                )

            # Also feed the FFT golden model with the new DDS sample
            self._fft_gold.push_dds(item.dds_amplitude)

            # When the FFT model has a full frame, compute the golden spectrum
            if self._fft_gold.ready():
                self._fft_golden_frame = self._fft_gold.compute()
                self._fft_out_idx      = 0
                self.logger.info(
                    "FFT golden model: full chirp frame received – "
                    "golden spectrum computed."
                )

        # ═══════════════════════════════════════════════════════════════════
        # STAGE 2 – FFT output vs. golden spectrum
        # ═══════════════════════════════════════════════════════════════════
        if item.fft_valid:
            if self._fft_out_idx < len(self._fft_golden_frame):
                ref_re, ref_im = self._fft_golden_frame[self._fft_out_idx]
                dut_re = _sign16(item.fft_re)
                dut_im = _sign16(item.fft_im)

                self._compare("FFT_re", dut_re, ref_re, TOLERANCE_FFT)
                self._compare("FFT_im", dut_im, ref_im, TOLERANCE_FFT)
                self._fft_out_idx += 1

                # When we have seen a full FFT output frame, load it into
                # the bit-reversal golden model as the reference frame.
                if self._fft_out_idx == N:
                    self._bitrev_gold.load_fft_frame(self._fft_golden_frame)
                    self._bitrev_out_idx = 0
            else:
                self.logger.warning(
                    f"FFT valid_out but no golden frame ready "
                    f"(fft_out_idx={self._fft_out_idx}) – skipping comparison."
                )

        # ═══════════════════════════════════════════════════════════════════
        # STAGE 3 – Bit-reversal output
        # ═══════════════════════════════════════════════════════════════════
        if item.bit_rev_valid:
            expected = self._bitrev_gold.get_expected(self._bitrev_out_idx)
            dut_re   = _sign16(item.bit_rev_re)
            dut_im   = _sign16(item.bit_rev_im)

            if expected is not None:
                ref_re, ref_im = expected
                self._compare("BITREV_re", dut_re, ref_re, TOLERANCE_BITREV)
                self._compare("BITREV_im", dut_im, ref_im, TOLERANCE_BITREV)
            else:
                self.logger.warning(
                    f"bit_rev_valid asserted but no reference ready "
                    f"(output_idx={self._bitrev_out_idx}) – skipping."
                )

            # Collect the bit-reversal frame so the MUX radar reference
            # can be built from the >>>7-scaled version of it.
            self._bitrev_frame.append((dut_re, dut_im))
            self._bitrev_out_idx += 1

            # After a full bit-reversal frame, build the MUX radar reference:
            # the RTL applies  >>>7  before writing into the radar RAM,
            # so our reference is  bit_rev_sample >> 7  (arithmetic shift).
            if self._bitrev_out_idx == N:
                # First 1667 samples of the bit-rev output go into the radar RAM
                radar_pairs = [
                    (re >> 7, im >> 7)
                    for (re, im) in self._bitrev_frame[:1667]
                ]
                self._mux_gold.load_radar(radar_pairs)
                self._bitrev_frame.clear()
                self._bitrev_out_idx = 0
                self.logger.info(
                    "Bit-reversal frame complete – radar reference loaded into MUX golden model."
                )

        # ═══════════════════════════════════════════════════════════════════
        # STAGE 4 – MUX output (frame structure + data)
        # ═══════════════════════════════════════════════════════════════════
        if item.mux_valid:
            idx      = self._mux_out_idx
            expected = self._mux_gold.get_expected(idx)
            dut_re   = _sign16(item.mux_re)
            dut_im   = _sign16(item.mux_im)

            if expected is not None:
                ref_re, ref_im, seg = expected
                self._compare("MUX_re", dut_re, ref_re, TOLERANCE_MUX,
                              extra=f"[{seg} bin {idx}]")
                self._compare("MUX_im", dut_im, ref_im, TOLERANCE_MUX,
                              extra=f"[{seg} bin {idx}]")
            else:
                self.logger.warning(
                    f"mux_valid asserted at sample_idx={idx} but no reference "
                    f"ready – skipping comparison."
                )

            # Feed MUX output into the TX IFFT golden model
            self._tx_gold.push_mux(item.mux_re, item.mux_im)

            self._mux_out_idx += 1
            if self._mux_out_idx == N:
                self._mux_out_idx = 0   # reset for next frame
                if self._tx_gold.ready():
                    self._tx_gold.compute()
                    self.logger.info(
                        "TX golden model: full MUX frame received – "
                        "golden TX output computed."
                    )

        # ═══════════════════════════════════════════════════════════════════
        # STAGE 5 – TX final output (tx_valid)
        # ═══════════════════════════════════════════════════════════════════
        if item.tx_valid:
            golden_pair = self._tx_gold.pop()

            if golden_pair is None:
                self.logger.warning(
                    "DUT asserted tx_valid but golden TX model has no output ready. "
                    "Possible frame-alignment issue – skipping comparison."
                )
                return

            ref_re, ref_im = golden_pair
            dut_re = _sign16(item.tx_out_re)
            dut_im = _sign16(item.tx_out_im)

            self._compare("TX_re", dut_re, ref_re, TOLERANCE_TX)
            self._compare("TX_im", dut_im, ref_im, TOLERANCE_TX)

            # Accumulate SQNR statistics
            self._sqnr_sig_pow   += ref_re**2 + ref_im**2
            self._sqnr_noise_pow += (dut_re - ref_re)**2 + (dut_im - ref_im)**2

    # ── _compare ──────────────────────────────────────────────────────────────
    def _compare(self, key: str, dut_val: int, ref_val: int,
                 tolerance: int, extra: str = ""):
        """Compare one signal with ±tolerance LSB window and update counters."""
        diff = abs(dut_val - ref_val)
        if diff <= tolerance:
            self.correct[key] += 1
        else:
            self.errors[key] += 1
            label = f"{key} {extra}".strip()
            self.logger.error(
                f"MISMATCH {label} "
                f"(error #{self.errors[key]}) "
                f"| DUT={dut_val}  REF={ref_val}  diff={diff}  tolerance={tolerance}"
            )

    # ── report_phase ──────────────────────────────────────────────────────────
    def report_phase(self):
        # SQNR
        if self._sqnr_noise_pow > 0:
            sqnr_db = 10 * np.log10(self._sqnr_sig_pow / self._sqnr_noise_pow)
        else:
            sqnr_db = float("inf")

        sqnr_ok  = sqnr_db >= SQNR_THRESHOLD_DB
        all_pass = all(v == 0 for v in self.errors.values()) and sqnr_ok
        overall  = "PASS ✓" if all_pass else "FAIL ✗"

        self.logger.info("╔══════════════════════════════════════════════════════╗")
        self.logger.info("║             TX PATH SCOREBOARD REPORT                ║")
        self.logger.info("╠══════════════════════════════════════════════════════╣")
        self.logger.info(f"║  Reset cycles detected        : {self.reset_cycles:<22} ║")
        self.logger.info("║ ─────────────────────────────────────────────────── ║")

        stage_rows = [
            ("DDS amplitude",   "DDS",       "DDS"),
            ("FFT real",        "FFT_re",    "FFT_re"),
            ("FFT imag",        "FFT_im",    "FFT_im"),
            ("Bit-rev real",    "BITREV_re", "BITREV_re"),
            ("Bit-rev imag",    "BITREV_im", "BITREV_im"),
            ("MUX real",        "MUX_re",    "MUX_re"),
            ("MUX imag",        "MUX_im",    "MUX_im"),
            ("TX out real",     "TX_re",     "TX_re"),
            ("TX out imag",     "TX_im",     "TX_im"),
        ]

        for label, c_key, e_key in stage_rows:
            total  = self.correct[c_key] + self.errors[e_key]
            self.logger.info(
                f"║  {label:<28}: {self.correct[c_key]}/{total:<16} errors={self.errors[e_key]:<5} ║"
            )

        self.logger.info("║ ─────────────────────────────────────────────────── ║")
        sqnr_str = f"{sqnr_db:.2f} dB" if sqnr_db != float('inf') else "∞ (perfect)"
        self.logger.info(f"║  TX output SQNR               : {sqnr_str:<22} ║")
        self.logger.info(f"║  SQNR threshold               : {SQNR_THRESHOLD_DB:.1f} dB{'':<17} ║")
        self.logger.info("║ ─────────────────────────────────────────────────── ║")
        self.logger.info(f"║  OVERALL RESULT               : {overall:<22} ║")
        self.logger.info("╚══════════════════════════════════════════════════════╝")

        if not all_pass:
            err_summary = ", ".join(
                f"{k}={v}" for k, v in self.errors.items() if v > 0
            )
            if not sqnr_ok:
                err_summary += f", SQNR={sqnr_db:.2f}dB<{SQNR_THRESHOLD_DB}dB"
            self.logger.critical(
                f"Scoreboard detected errors: {err_summary}"
            )