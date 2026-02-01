onerror {resume}
quietly WaveActivateNextPane {} 0
add wave -noupdate -expand -group DUT_signals -radix unsigned /ALU_top/DUT/a
add wave -noupdate -expand -group DUT_signals -radix unsigned /ALU_top/DUT/b
add wave -noupdate -expand -group DUT_signals -radix unsigned /ALU_top/DUT/op
add wave -noupdate -expand -group DUT_signals -radix unsigned /ALU_top/DUT/out
add wave -noupdate -expand -group DUT_signals /ALU_top/DUT/c
add wave -noupdate -expand -group Reference_signals -color Orange -radix unsigned /ALU_scoreboard_pkg::out_ref
add wave -noupdate -expand -group Reference_signals -color Orange -radix unsigned /ALU_scoreboard_pkg::c_ref
TreeUpdate [SetDefaultTree]
WaveRestoreCursors {{Cursor 1} {0 ns} 0}
quietly wave cursor active 1
configure wave -namecolwidth 226
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
WaveRestoreZoom {0 ns} {33 ns}
