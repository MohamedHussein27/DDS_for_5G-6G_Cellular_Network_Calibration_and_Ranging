# SYSTEM_TOP UVM Verification Environment

**Sponsor:** Analog Devices, Inc. (ADI)
**Institution:** Faculty of Engineering, Ain Shams University
**Project:** DDS for 5G/6G Cellular Network Calibration and Ranging

---

## Overview

Full UVM-pyuvm environment for the `SYSTEM_TOP` ISAC DUT — the
transmitter (TX_TOP) wired directly to the receiver (RX_TOP) with no
channel in between.  Three output paths are verified in a single
simulation:

| Path  | DUT signals                               | Frame size |
|-------|-------------------------------------------|------------|
| TX    | `tx_valid`, `tx_out_re`, `tx_out_im`      | 4 096      |
| OFDM  | `ofdm_valid_out`, `ofdm_out_re/im`        | 2 048      |
| Radar | `radar_valid_out`, `radar_out_re/im`      | 2 048      |

---

## File Map

```
system_top_env/
├── system_top_seq_item.py      Unified TX+RX transaction packet
├── system_top_sequencer.py     Pass-through sequencer
├── system_top_driver.py        Drives TX bus + backdoor ROM writes
├── system_top_monitor.py       Samples all DUT boundary signals
├── system_top_agent.py         Active/passive top-level agent
├── tx_agent.py                 Active/passive TX sub-agent (re-uses TX standalone)
├── rx_agent.py                 Active/passive RX sub-agent (re-uses RX standalone)
├── system_top_subscriber.py    Functional coverage (21 cover points/crosses)
├── system_top_scoreboard.py    Three-path golden-model scoreboard (provided)
├── system_top_golden_model.py  Full TX→RX system golden pipeline (provided)
├── system_top_env.py           Environment: wires agents→scoreboard→subscriber
├── system_top_test.py          Base test + chirp test + random test + sequences
└── Makefile                    cocotb-make build
```

---

## Architecture

```
                           ┌─────────────────────────────────┐
                           │         system_top_env           │
                           │                                   │
  ┌─────────────────────┐  │  ┌─────────────────────────────┐ │
  │  system_top_agent   │  │  │  system_top_scoreboard      │ │
  │  (active/passive)   │──┼──│  TX + OFDM + Radar paths    │ │
  │  sqr → drv → mon   │  │  └─────────────────────────────┘ │
  │  agt_ap ────────────┼──┼──► subscriber (coverage)         │
  └─────────────────────┘  │                                   │
                           │  ┌───────────────┐               │
  ┌─────────────────────┐  │  │  tx_agent     │               │
  │  tx_agent           │──┼──► dds_export    │               │
  │  (active/passive)   │  │  └───────────────┘               │
  │  top_driver/monitor │  │                                   │
  └─────────────────────┘  │  ┌───────────────┐               │
                           │  │  rx_agent     │               │
  ┌─────────────────────┐  │  │  (passive by  │               │
  │  rx_agent           │──┼──► default)      │               │
  │  (active/passive)   │  │  └───────────────┘               │
  │  rx_driver/monitor  │  │                                   │
  └─────────────────────┘  └─────────────────────────────────┘
```

---

## Agent Activity Modes

All three agents are independently configurable via `ConfigDB`:

```python
# Default for loopback (SYSTEM_TOP with internal TX→RX wire)
ConfigDB().set(None, "uvm_test_top.env.sys_agent.*", "is_active", UVM_ACTIVE)
ConfigDB().set(None, "uvm_test_top.env.tx_agent.*",  "is_active", UVM_ACTIVE)
ConfigDB().set(None, "uvm_test_top.env.rx_agent.*",  "is_active", UVM_PASSIVE)

# To inject independent RX stimulus (e.g. channel injection test):
ConfigDB().set(None, "uvm_test_top.env.rx_agent.*",  "is_active", UVM_ACTIVE)
```

---

## Running Tests

```bash
# Directed chirp test (default)
make chirp

# Constrained-random test
make random

# Explicit test name
make TEST=system_top_chirp_test

# Clean sim artefacts
make clean
```

---

## Sequence Library (`system_top_test.py`)

| Sequence            | Purpose                                               |
|---------------------|-------------------------------------------------------|
| `reset_seq`         | Assert `rst_n=0` for `RESET_CYCLES`, then de-assert  |
| `backdoor_ofdm_seq` | Zero-time write to internal OFDM ROM                  |
| `program_chirp_seq` | Write FTW registers + trigger `dds_enable`            |
| `idle_seq`          | Hold bus idle for N cycles while DUT computes         |

---

## Scoreboard Logic

The `system_top_scoreboard` calls `run_system_top_pipeline()` **once**
on the first `dds_item.enable` pulse to produce all six reference
vectors.  No sample accumulation is needed because the golden model
only requires `(FTW_start, FTW_step, cycles, ofdm_re[], ofdm_im[])`.

The OFDM ROM contents are read via DUT backdoor handle exactly as in
the standalone TX scoreboard.

---

## Functional Coverage

21 cover points and crosses in `system_top_subscriber.py`:

- Control signals: `rst_n`, `wr_en`, `rd_en`, `addr`, `dds_ready_flag`,
  `tx_valid`, `rx_valid_in`, `ofdm_valid_out`, `radar_valid_out`
- Data ranges (min_neg / negative / zero / positive / max_pos) for all
  six 16-bit output signals
- System-level cross: `tx_valid × ofdm_valid_out × radar_valid_out`
  (ensures all three output paths fire simultaneously)

Coverage data exported to `system_top_coverage.xml` in `report_phase`.

---

## Python Path Dependencies

These files must be on `PYTHONPATH` at runtime (already handled by the
Makefile):

```
# From TX standalone env:
top_seq_item.py, top_driver.py, top_monitor.py, top_sequencer.py

# From RX standalone env:
rx_top_driver.py, rx_top_monitor.py

# Golden model chain:
dds_new_model.py, fft_fixed.py, ifft_fixed.py, mux_new.py, rx_golden.py
```
