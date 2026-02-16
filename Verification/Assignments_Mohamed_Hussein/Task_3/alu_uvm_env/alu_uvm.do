vlib work

vlog alu_shared_pkg.sv alu_config.sv alu_interface.sv dut.sv \
     alu_sequence_item.sv alu_sequencer.sv alu_driver.sv \
     alu_monitor.sv alu_agent.sv alu_scoreboard.sv \
     alu_subscriber.sv alu_env.sv \
     alu_reset_sequence.sv alu_main_sequence.sv \
     alu_add_sequence.sv alu_xor_sequence.sv \
     alu_and_sequence.sv alu_or_sequence.sv \
     alu_test.sv alu_top.sv \
     +cover -covercells

########################################
# TEST 1
########################################
vsim -voptargs=+acc work.alu_top -cover +UVM_TESTNAME=alu_add_xor_test

run 0

add wave -position insertpoint \
sim:/alu_top/dut/clk \
sim:/alu_top/dut/rst_n \
sim:/alu_top/dut/a \
sim:/alu_top/dut/b \
sim:/alu_top/dut/op \
sim:/alu_top/dut/out \
sim:/alu_top/dut/c \
sim:/alu_shared_pkg::error_count_out \
sim:/alu_shared_pkg::error_count_c
run -all
coverage save alu_add_xor.ucdb
quit -sim


########################################
# TEST 2
########################################
vsim -voptargs=+acc work.alu_top -cover +UVM_TESTNAME=alu_and_or_test

run 0

add wave -position insertpoint \
sim:/alu_top/dut/clk \
sim:/alu_top/dut/rst_n \
sim:/alu_top/dut/a \
sim:/alu_top/dut/b \
sim:/alu_top/dut/op \
sim:/alu_top/dut/out \
sim:/alu_top/dut/c \
sim:/alu_shared_pkg::error_count_out \
sim:/alu_shared_pkg::error_count_c
run -all
coverage save alu_and_or.ucdb
quit -sim


########################################
# TEST 3
########################################
vsim -voptargs=+acc work.alu_top -cover +UVM_TESTNAME=alu_rand_test

run 0

add wave -position insertpoint \
sim:/alu_top/dut/clk \
sim:/alu_top/dut/rst_n \
sim:/alu_top/dut/a \
sim:/alu_top/dut/b \
sim:/alu_top/dut/op \
sim:/alu_top/dut/out \
sim:/alu_top/dut/c \
sim:/alu_shared_pkg::error_count_out \
sim:/alu_shared_pkg::error_count_c

run -all
coverage save alu_rand.ucdb

########################################
# MERGE COVERAGE
########################################
vcover merge alu_merged.ucdb alu_add_xor.ucdb alu_and_or.ucdb alu_rand.ucdb

vcover report alu_merged.ucdb
