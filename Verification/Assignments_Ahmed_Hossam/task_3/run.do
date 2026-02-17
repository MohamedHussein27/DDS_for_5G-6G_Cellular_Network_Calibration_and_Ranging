vlib work
vlog -f sourcefile.txt +cover -covercells
vsim -voptargs=+acc work.top -cover +UVM_TESTNAME=test_all
add wave -r /*
coverage save alu_tb.ucdb -onexit
run -all
