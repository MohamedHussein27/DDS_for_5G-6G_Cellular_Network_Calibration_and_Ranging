"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_sequences.py

    Description:
        This module contains a comprehensive suite of UVM stimulus sequences 
        for verifying the 4096-point IFFT block. It translates mathematical 
        frequency-domain matrices into cycle-accurate, 16-bit signed fixed-point 
        (Q11.5 format) transactions. The suite includes corner cases (max positive/
        negative, resets), protocol checks (bursts, back-to-back frames), and 
        fundamental DSP signal models (impulses, single-tone, and two-tone 
        generation) for exhaustive functional coverage.
"""

import cocotb
from cocotb.triggers import *
from pyuvm import *
import logging
import numpy as np
from ifft_item import ifft_item

# ─────────────────────────────────────────────────────────────────────────────
# RTL / Q8.8 constants  (shared across all sequences)
# ─────────────────────────────────────────────────────────────────────────────
N       = 4096
WL      = 16
Q_FRAC  = 5
MAX_VAL = (1 << (WL - 1)) - 1   #  32767
MIN_VAL = -(1 << (WL - 1))      # -32768

def _f2q(f):
    """float → signed Q8.8 integer, saturated to [MIN_VAL, MAX_VAL]."""
    return max(MIN_VAL, min(MAX_VAL, int(round(f * (1 << Q_FRAC)))))

def _frame(arr):
    """numpy complex array → list of (re_q8, im_q8) signed int pairs."""
    return [(_f2q(c.real), _f2q(c.imag)) for c in arr]


# ─────────────────────────────────────────────────────────────────────────────
# 1. reset_sequence
# ─────────────────────────────────────────────────────────────────────────────
class reset_sequence(uvm_sequence):
    def __init__(self, name="reset_sequence"):
        super().__init__(name)

    async def body(self):
        self.seq_item = ifft_item.create("seq_item")
        await self.start_item(self.seq_item)
        self.seq_item.randomize_with(lambda rst_n: rst_n in [0])
        await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 2. allzero_sequence  – all bins 0+0j  →  all-zero output
# ─────────────────────────────────────────────────────────────────────────────
class allzero_sequence(uvm_sequence):
    def __init__(self, name="allzero_sequence"):
        super().__init__(name)

    async def body(self):
        for _ in range(N):
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag:
                    rst_n == 1 and valid_in == 1 and in_real == 0 and in_imag == 0
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 3. impulse_sequence  – X[0] = 1+0j, rest = 0  →  constant 1/N output
# ─────────────────────────────────────────────────────────────────────────────
class impulse_sequence(uvm_sequence):
    def __init__(self, name="impulse_sequence"):
        super().__init__(name)

    async def body(self):
        X = np.zeros(N, dtype=complex)
        X[0] = 1.0 + 0j
        pairs = _frame(X)

        for re, im in pairs:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag, _re=re, _im=im:
                    rst_n == 1 and valid_in == 1 and in_real == _re and in_imag == _im
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 4. dc_sequence  – all bins = 1+0j  →  large spike at n=0, zeros elsewhere
# ─────────────────────────────────────────────────────────────────────────────
class dc_sequence(uvm_sequence):
    def __init__(self, name="dc_sequence"):
        super().__init__(name)

    async def body(self):
        re_val = _f2q(1.0)
        im_val = 0

        for _ in range(N):
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag, _re=re_val, _im=im_val:
                    rst_n == 1 and valid_in == 1 and in_real == _re and in_imag == _im
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 5. singletone_sequence  – energy at one bin k  →  pure sinusoid at k/N
# ─────────────────────────────────────────────────────────────────────────────
class singletone_sequence(uvm_sequence):
    def __init__(self, name="singletone_sequence", k=1, amplitude=1.0):
        super().__init__(name)
        self.k         = k
        self.amplitude = amplitude

    async def body(self):
        X = np.zeros(N, dtype=complex)
        X[self.k] = self.amplitude + 0j
        pairs = _frame(X)

        for re, im in pairs:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag, _re=re, _im=im:
                    rst_n == 1 and valid_in == 1 and in_real == _re and in_imag == _im
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 6. twotone_sequence  – bins k1 and k2  →  superposition of two sinusoids
# ─────────────────────────────────────────────────────────────────────────────
class twotone_sequence(uvm_sequence):
    def __init__(self, name="twotone_sequence", k1=10, k2=200, a1=0.5, a2=0.5):
        super().__init__(name)
        self.k1, self.k2 = k1, k2
        self.a1, self.a2 = a1, a2

    async def body(self):
        X = np.zeros(N, dtype=complex)
        X[self.k1] = self.a1
        X[self.k2] = self.a2
        pairs = _frame(X)

        for re, im in pairs:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag, _re=re, _im=im:
                    rst_n == 1 and valid_in == 1 and in_real == _re and in_imag == _im
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 7. maxpositive_sequence  – every sample in_real = +32767
# ─────────────────────────────────────────────────────────────────────────────
class maxpositive_sequence(uvm_sequence):
    def __init__(self, name="maxpositive_sequence"):
        super().__init__(name)

    async def body(self):
        for _ in range(N):
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag:
                    rst_n == 1 and valid_in == 1 and in_real == MAX_VAL and in_imag == 0
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 8. maxnegative_sequence  – every sample in_real = -32768
# ─────────────────────────────────────────────────────────────────────────────
class maxnegative_sequence(uvm_sequence):
    def __init__(self, name="maxnegative_sequence"):
        super().__init__(name)

    async def body(self):
        for _ in range(N):
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag:
                    rst_n == 1 and valid_in == 1 and in_real == MIN_VAL and in_imag == 0
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 9. alternating_sequence  – +MAX / -MAX every other sample
# ─────────────────────────────────────────────────────────────────────────────
class alternating_sequence(uvm_sequence):
    def __init__(self, name="alternating_sequence"):
        super().__init__(name)

    async def body(self):
        for i in range(N):
            val = MAX_VAL if i % 2 == 0 else MIN_VAL
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag, _v=val:
                    rst_n == 1 and valid_in == 1 and in_real == _v and in_imag == 0
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 10. halfnyquist_sequence  – bin N/2 = 2048  →  alternating ±1/N output
# ─────────────────────────────────────────────────────────────────────────────
class halfnyquist_sequence(uvm_sequence):
    def __init__(self, name="halfnyquist_sequence"):
        super().__init__(name)

    async def body(self):
        X = np.zeros(N, dtype=complex)
        X[N // 2] = 1.0 + 0j
        pairs = _frame(X)

        for re, im in pairs:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag, _re=re, _im=im:
                    rst_n == 1 and valid_in == 1 and in_real == _re and in_imag == _im
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 11. ramp_sequence  – in_real ramps 0 → ~0.9 in Q8.8, in_imag = 0
# ─────────────────────────────────────────────────────────────────────────────
class ramp_sequence(uvm_sequence):
    def __init__(self, name="ramp_sequence", max_float=0.9):
        super().__init__(name)
        self.max_float = max_float

    async def body(self):
        ramp = [_f2q(f) for f in np.linspace(0.0, self.max_float, N)]

        for val in ramp:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag, _v=val:
                    rst_n == 1 and valid_in == 1 and in_real == _v and in_imag == 0
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 12. random_sequence  – fully random in_real and in_imag, rst_n forced to 1
# ─────────────────────────────────────────────────────────────────────────────
class random_sequence(uvm_sequence):
    def __init__(self, name="random_sequence"):
        super().__init__(name)

    async def body(self):
        for _ in range(N):
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in: rst_n == 1 and valid_in == 1
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 13. burst_sequence  – valid_in=1 for burst_len samples, then idle
# ─────────────────────────────────────────────────────────────────────────────
class burst_sequence(uvm_sequence):
    def __init__(self, name="burst_sequence", burst_len=512):
        super().__init__(name)
        self.burst_len = burst_len

    async def body(self):
        # Active burst — random data, valid_in=1
        for _ in range(self.burst_len):
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in: rst_n == 1 and valid_in == 1
            )
            await self.finish_item(self.seq_item)

        # Idle for the rest of the frame window — valid_in=0, inputs=0
        for _ in range(N - self.burst_len):
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag:
                    rst_n == 1 and valid_in == 0 and in_real == 0 and in_imag == 0
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 14. resetmidframe_sequence  – rst_n=0 injected mid-frame, then clean frame
# ─────────────────────────────────────────────────────────────────────────────
class resetmidframe_sequence(uvm_sequence):
    def __init__(self, name="resetmidframe_sequence", reset_at=N // 2, reset_cycles=5):
        super().__init__(name)
        self.reset_at     = reset_at
        self.reset_cycles = reset_cycles

    async def body(self):
        # Dirty half-frame (single tone at bin 1)
        X_dirty = np.zeros(N, dtype=complex)
        X_dirty[1] = 0.3 + 0j
        pairs_dirty = _frame(X_dirty)

        for re, im in pairs_dirty[:self.reset_at]:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag, _re=re, _im=im:
                    rst_n == 1 and valid_in == 1 and in_real == _re and in_imag == _im
            )
            await self.finish_item(self.seq_item)

        # Assert reset for reset_cycles cycles
        for _ in range(self.reset_cycles):
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(lambda rst_n: rst_n in [0])
            await self.finish_item(self.seq_item)

        # One idle after reset release so RTL stabilises
        self.seq_item = ifft_item.create("seq_item")
        await self.start_item(self.seq_item)
        self.seq_item.randomize_with(
            lambda rst_n, valid_in, in_real, in_imag:
                rst_n == 1 and valid_in == 0 and in_real == 0 and in_imag == 0
        )
        await self.finish_item(self.seq_item)

        # Clean full frame — impulse at bin 0 (easy golden model check)
        X_clean = np.zeros(N, dtype=complex)
        X_clean[0] = 1.0 + 0j
        pairs_clean = _frame(X_clean)

        for re, im in pairs_clean:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag, _re=re, _im=im:
                    rst_n == 1 and valid_in == 1 and in_real == _re and in_imag == _im
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 15. backtoback_sequence  – two consecutive frames, valid_in stays high
# ─────────────────────────────────────────────────────────────────────────────
class backtoback_sequence(uvm_sequence):
    def __init__(self, name="backtoback_sequence", seed=7):
        super().__init__(name)
        self.seed = seed

    async def body(self):
        # Frame 1: single tone at bin 5
        X1 = np.zeros(N, dtype=complex)
        X1[5] = 0.4 + 0j
        pairs1 = _frame(X1)

        # Frame 2: pre-computed random values
        rng    = np.random.default_rng(self.seed)
        re_v   = rng.uniform(-0.3, 0.3, N)
        im_v   = rng.uniform(-0.3, 0.3, N)
        pairs2 = [(_f2q(r), _f2q(i)) for r, i in zip(re_v, im_v)]

        for re, im in pairs1 + pairs2:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in, in_real, in_imag, _re=re, _im=im:
                    rst_n == 1 and valid_in == 1 and in_real == _re and in_imag == _im
            )
            await self.finish_item(self.seq_item)


# ─────────────────────────────────────────────────────────────────────────────
# 16. fullregression_sequence  – runs every sequence in order
# ─────────────────────────────────────────────────────────────────────────────
class fullregression_sequence(uvm_sequence):
    def __init__(self, name="fullregression_sequence"):
        super().__init__(name)

    async def body(self):
        sequences = [
            reset_sequence        ("reset"),
            allzero_sequence      ("allzero"),
            impulse_sequence      ("impulse"),
            dc_sequence           ("dc"),
            singletone_sequence   ("tone_k1",    k=1),
            singletone_sequence   ("tone_k100",  k=100),
            singletone_sequence   ("tone_k2047", k=2047),
            twotone_sequence      ("twotone"),
            maxpositive_sequence  ("maxpos"),
            maxnegative_sequence  ("maxneg"),
            alternating_sequence  ("alternating"),
            halfnyquist_sequence  ("nyquist"),
            ramp_sequence         ("ramp"),
            random_sequence       ("random"),
            burst_sequence        ("burst_512",  burst_len=512),
            resetmidframe_sequence("rst_mid"),
            backtoback_sequence   ("b2b"),
        ]

        for seq in sequences:
            self.logger.info(f"[fullregression_sequence] ── Starting {seq.get_name()} ──")
            await seq.start(self.sequencer)