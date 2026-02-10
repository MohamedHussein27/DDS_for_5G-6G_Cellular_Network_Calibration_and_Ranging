onerror {resume}
quietly WaveActivateNextPane {} 0
add wave -noupdate -expand -group {DUT signals} -radix unsigned /ALU_top/dut/clk
add wave -noupdate -expand -group {DUT signals} -radix unsigned /ALU_top/dut/a
add wave -noupdate -expand -group {DUT signals} -radix unsigned /ALU_top/dut/b
add wave -noupdate -expand -group {DUT signals} -radix unsigned /ALU_top/dut/op
add wave -noupdate -expand -group {DUT signals} -radix unsigned /ALU_top/dut/c
add wave -noupdate -expand -group {DUT signals} -radix unsigned /ALU_top/dut/out
add wave -noupdate -expand -group {DUT signals} -radix unsigned /ALU_top/dut/out_data
TreeUpdate [SetDefaultTree]
WaveRestoreCursors {{Cursor 1} {2 ns} 0}
quietly wave cursor active 1
configure wave -namecolwidth 150
configure wave -valuecolwidth 100
configure wave -justifyvalue left
configure wave -signalnamewidth 1
configure wave -snapdistance 10
configure wave -datasetprefix 0
configure wave -rowmargin 4
configure wave -childrowmargin 2
configure wave -gridoffset 0
configure wave -gridperiod 1
configure wave -griddelta 40
configure wave -timeline 0
configure wave -timelineunits ns
update
WaveRestoreZoom {2082 ns} {2101 ns}
