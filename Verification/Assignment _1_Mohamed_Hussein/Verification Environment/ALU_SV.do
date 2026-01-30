vlib work
vlog DUT.sv ALU_IF.sv ALU_SHARED.sv ALU_SEQUENCE_ITEM.sv ALU_COVERAGE.sv ALU_MONITOR.sv ALU_SCOREBOARD.sv ALU_DRIVER.sv ALU_TOP.sv +cover -covercells
vsim -voptargs=+acc work.alu_top -cover 
add wave -position insertpoint  \
sim:/alu_top/dut/a \
sim:/alu_top/dut/b \
sim:/alu_top/dut/op \
sim:/alu_top/dut/c \
sim:/alu_top/dut/out \
sim:/alu_top/dut/out_data \

add wave -position insertpoint  \
sim:/alu_top/MONITOR/tr_mon.a \
sim:/alu_top/MONITOR/tr_mon.b \
sim:/alu_top/MONITOR/tr_mon.op \
sim:/alu_top/MONITOR/tr_mon.out \
sim:/alu_top/MONITOR/tr_mon.c \

add wave -position insertpoint  \
sim:/alu_top/MONITOR/scb_mon.out_ref \
sim:/alu_top/MONITOR/scb_mon.c_ref

add wave -position insertpoint  \
sim:/shared_pkg::error_count_out \
sim:/shared_pkg::correct_count_out \
sim:/shared_pkg::error_count_c \
sim:/shared_pkg::correct_count_c
coverage save ALU.ucdb -onexit 
run -all