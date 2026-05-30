# ==============================================================================
# cocotb Makefile: ADI 5G/6G Cellular Network Calibration - RX Path (ISAC)
# ==============================================================================
export COCOTB_RESOLVE_X = ZEROS


# Simulator Defaults
GUI ?= 1
WAVES ?= 1
SIM ?= questa
TOPLEVEL_LANG ?= verilog



# VERILOG SOURCES 
# ─────────────────────────────────────────────────────────
V_FILES = bit_reversal.v \
          butterfly.v \
          config_reg.v \
          controlunit_2048.v \
          controlunit_2048pts.v \
          controlunit_4096.v \
          controlunit_4096pts.v \
          DDS.v \
          delayfeedback.v \
          fft__4096_stage.v \
          fft_4096_stage.v \
          fft_4096_top.v \
          first_quad_address.v \
          ifft_2048_stage.v \
          ifft_2048_top.v \
          ifft_4096_stage.v \
          ifft_4096_top.v \
          isac_demux.v \
          isac_mux.v \
          LUT.v \
          multiplier.v \
          multiplier_RX.v \
          negative_mux.v \
          ofdm_rom.v \
          PHASE_ACC.v \
          RAM_RX.v \
          RX_top.v \
          SYSTEM_TOP.v \
          tb_fft_4096.v \
          twiddlerom_2048.v \
          twiddlerom_4096.v \
          twiddlerom_4096pts.v \
          TX_TOP.v

# $(CURDIR) ensures the simulator finds the files from any directory
VERILOG_SOURCES = $(foreach f, $(V_FILES), $(CURDIR)/$(f))

# ─────────────────────────────────────────────────────────
# Top-level Verilog module name
TOPLEVEL = SYSTEM_TOP

# Python test module name (Default fallback, overridden by Makefile.runs)
MODULE = rx_top_test

# Tell the Questa compiler to instrument the RTL for coverage (bcest)
COMPILE_ARGS += +cover=bcest

include $(shell cocotb-config --makefiles)/Makefile.sim

# ─────────────────────────────────────────────────────────
# CLEAN TARGET (Must use a real TAB on the line below!)
clean::
	rm -rf __pycache__ results.xml modelsim.ini transcript vsim.wlf sim_build

