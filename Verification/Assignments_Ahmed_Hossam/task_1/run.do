vlib work
vdel -all
vlog -sv +cover alu_pack.sv
vlog -sv +cover alu_if.sv
vlog -sv +cover DUT.sv
vlog -sv alu_env.sv
vlog -sv tb_top.sv

vsim -voptargs=+acc work.tb_top -cover
add wave -r /*
coverage save alu_tb.ucdb -onexit
run -all
