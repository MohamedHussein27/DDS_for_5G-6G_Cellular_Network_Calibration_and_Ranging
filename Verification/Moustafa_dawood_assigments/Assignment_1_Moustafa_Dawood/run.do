vlib work
vdel -all

vlog -sv +cover alu_pkg.sv
vlog -sv +cover alu_if.sv
vlog -sv +cover alu.sv      ;# DUT
vlog -sv alu_env.sv
vlog -sv alu_top.sv

vsim -coverage top
run -all

coverage report -summary
coverage report -codeAll -details -output coverage_report.txt
