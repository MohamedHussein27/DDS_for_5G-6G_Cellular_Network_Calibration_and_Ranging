# ─────────────────────────────────────────────────────────────────────────────
# 15. backtoback_sequence  
# 
# WHAT IT DOES:
# =============
# This sequence drives two consecutive frames back-to-back WITHOUT any idle
# cycles in between (valid_in stays high continuously). The frames are:
#
#   FRAME 1: Single tone at bin 5 (amplitude 0.4) → produces a clean sinusoid
#   FRAME 2: Random values (uniform between -0.3 and 0.3) → stresses the core
#            with unpredictable data
#
# The key feature is that pairs1 + pairs2 concatenates both frames so the
# second frame starts immediately after the last sample of the first frame,
# with NO valid_in=0 gap.
#
# PURPOSE:
# ========
# Verify that:
#   - The IFFT core can handle continuous streaming without frame gaps
#   - Internal pipelines don't stall or corrupt data when back-to-back frames
#     are driven
#   - No loss of synchronization between frames
#   - The core correctly processes the first frame's output while simultaneously
#     receiving the second frame's input
#   - The random second frame doesn't cause overflow or unexpected behavior
# ─────────────────────────────────────────────────────────────────────────────
import numpy as np
from pyuvm import *
from ifft_item import ifft_item
from ifft_seq_utils import N, _frame, _f2q
class backtoback_sequence(uvm_sequence):
    def __init__(self, name="backtoback_sequence", seed=7):
        super().__init__(name)
        self.seed = seed    # Fixed seed for reproducible random values
        self.multi_frame = False
        self.latency = 12

    async def body(self):
        # ─────────────────────────────────────────────────────────────
        # FRAME 1: Single tone at bin 5 (amplitude 0.4)
        # Produces a clean, predictable sinusoid output
        # ─────────────────────────────────────────────────────────────
        X1 = np.zeros(N, dtype=complex)
        X1[5] = 0.4 + 0j
        pairs1 = _frame(X1)

        # ─────────────────────────────────────────────────────────────
        # FRAME 2: Random values with fixed seed (reproducible)
        # Values range from -0.3 to 0.3 to keep within reasonable bounds
        # ─────────────────────────────────────────────────────────────
        rng    = np.random.default_rng(self.seed)
        re_v   = rng.uniform(-0.3, 0.3, N)   # Random real parts
        im_v   = rng.uniform(-0.3, 0.3, N)   # Random imaginary parts
        pairs2 = [(_f2q(r), _f2q(i)) for r, i in zip(re_v, im_v)]

        # ─────────────────────────────────────────────────────────────
        # Back-to-back drive: FRAME 1 followed immediately by FRAME 2
        # NO idle cycles between frames (valid_in stays high throughout)
        # ─────────────────────────────────────────────────────────────
        for re, im in pairs1 + pairs2:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in:
                    rst_n == 1 and valid_in == 1
            )
            self.seq_item.in_real = re
            self.seq_item.in_imag = im
            await self.finish_item(self.seq_item)

        if not self.multi_frame:
            for _ in range(N + self.latency):
                self.seq_item = ifft_item.create("seq_item")
                await self.start_item(self.seq_item)
                self.seq_item.randomize_with(
                    lambda rst_n, valid_in: rst_n == 1 and valid_in == 0
                )
                self.seq_item.in_real = 0
                self.seq_item.in_imag = 0
                await self.finish_item(self.seq_item)