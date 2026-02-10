vlib work
vdel -all
vlog -sv  counter_if.sv
vlog -sv  count_fsm.svp
vlog -sv counter_env.sv
vlog -sv counter_top.sv

# Added -sv_seed random here
vsim -voptargs=+acc  work.tb_top -cover

add wave -r /*

run -all
coverage save counter_tb.ucdb -onexit
vcover report counter_tb.ucdb \
  -details \
  -annotate \
  -all \
  -output coverage.txt