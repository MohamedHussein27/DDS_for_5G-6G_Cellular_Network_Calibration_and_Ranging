vlib work

vlog count_fsm.svp COUNTER_IF.sv COUNTER_SHARED.sv COUNTER_SEQ_ITEM.sv COUNTER_COVERAGE.sv \
     COUNTER_MONITOR.sv COUNTER_SCOREBOARD.sv COUNTER_DRIVER.sv COUNTER_GENERATOR.sv COUNTER_AGENT.sv COUNTER_ENV.sv COUNTER_TOP.sv \
     +cover=sbt -covercells

vsim -voptargs=+acc work.counter_top -cover

add wave -position insertpoint  sim:/counter_top/dut/clk
add wave -position insertpoint  sim:/counter_top/dut/rst_n 
add wave -position insertpoint  sim:/counter_top/dut/start 
add wave -position insertpoint  sim:/counter_top/dut/wait_timer 
add wave -position insertpoint  sim:/counter_top/dut/flag 
add wave -position insertpoint  sim:/counter_top/dut/busy 
add wave -position insertpoint  sim:/counter_top/dut/count_value



add wave -position insertpoint  sim:/shared_pkg::cs 
add wave -position insertpoint  sim:/shared_pkg::ns 
add wave -position insertpoint  sim:/shared_pkg::fsm_out
add wave -position insertpoint  sim:/shared_pkg::deassert_count
add wave -position insertpoint  sim:/shared_pkg::internal_counter
add wave -position insertpoint  sim:/shared_pkg::count_ref  
add wave -position insertpoint  sim:/shared_pkg::error_count_out 
add wave -position insertpoint  sim:/shared_pkg::error_count_busy  


run -all
coverage save COUNTER.ucdb


vcover report COUNTER.ucdb \
  -details \
  -annotate \
  -all \
  -output coverage.txt
