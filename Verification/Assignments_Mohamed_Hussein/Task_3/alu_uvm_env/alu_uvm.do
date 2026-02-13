vlib work

vlog alu_shared_pkg.sv alu_config.sv alu_interface.sv dut.sv alu_sequence_item.sv alu_sequencer.sv alu_driver.sv alu_monitor.sv alu_agent.sv alu_scoreboard.sv alu_subscriber.sv alu_env.sv alu_reset_sequence.sv alu_main_sequence.sv alu_add_sequence.sv alu_xor_sequence.sv alu_and_sequence.sv alu_or_sequence.sv alu_test.sv alu_top.sv +cover -covercells

vsim -voptargs=+acc work.alu_top -cover

add wave -position insertpoint sim:/alu_top/dut/*

coverage save alu_test.ucdb -onexit

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