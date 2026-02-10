vlib work
vlog -f src_files.txt +cover -covercells
vsim -voptargs=+acc work.count_top -cover 
do wave.do
run -all
coverage save count_top.ucdb -onexit 
vcover report count_top.ucdb -details -annotate -all

