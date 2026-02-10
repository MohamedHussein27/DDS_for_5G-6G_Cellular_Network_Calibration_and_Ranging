vlib work
vlog -f src_files.txt +cover -covercells
vsim -voptargs=+acc work.ALU_top -cover 
do wave.do
run -all
coverage save ALU_top.ucdb -onexit 
vcover report ALU_top.ucdb -details -annotate -all

