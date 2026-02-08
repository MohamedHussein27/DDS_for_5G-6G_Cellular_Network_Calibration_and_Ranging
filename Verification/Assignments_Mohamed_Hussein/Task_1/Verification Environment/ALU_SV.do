vlib work

vlog DUT.sv ALU_IF.sv ALU_SHARED.sv ALU_SEQUENCE_ITEM.sv ALU_COVERAGE.sv \
     ALU_MONITOR.sv ALU_SCOREBOARD.sv ALU_DRIVER.sv ALU_GENERATOR.sv ALU_AGENT.sv ALU_ENV.sv ALU_TOP.sv \
     +cover=sbt -covercells

vsim -voptargs=+acc work.alu_top -cover

add wave -position insertpoint  sim:/alu_top/dut/clk
add wave -position insertpoint  sim:/alu_top/dut/rst_n 
add wave -position insertpoint  sim:/alu_top/dut/a 
add wave -position insertpoint  sim:/alu_top/dut/b 
add wave -position insertpoint  sim:/alu_top/dut/op 
add wave -position insertpoint  sim:/alu_top/dut/c 
add wave -position insertpoint  sim:/alu_top/dut/out 
add wave -position insertpoint  sim:/alu_top/dut/out_data 


add wave -position insertpoint  sim:/shared_pkg::error_count_out 
add wave -position insertpoint  sim:/shared_pkg::error_count_c 


run -all
coverage save ALU.ucdb


vcover report ALU.ucdb \
  -details \
  -annotate \
  -all \
  -output coverage.txt
