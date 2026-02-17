onerror {resume}
quietly WaveActivateNextPane {} 0
add wave -noupdate -expand -group DUT_signals -radix unsigned /count_top/dut/rst_n
add wave -noupdate -expand -group DUT_signals -radix unsigned /count_top/dut/clk
add wave -noupdate -expand -group DUT_signals -radix unsigned /count_top/dut/wait_timer
add wave -noupdate -expand -group DUT_signals -radix unsigned /count_top/dut/start
add wave -noupdate -expand -group DUT_signals -radix unsigned /count_top/dut/flag
add wave -noupdate -expand -group DUT_signals -radix unsigned /count_top/dut/count_value
add wave -noupdate -expand -group DUT_signals -radix unsigned /count_top/dut/busy
add wave -noupdate -expand -group ref_Signals -radix unsigned /count_scoreboard_pkg::next_state
add wave -noupdate -expand -group ref_Signals -radix unsigned /count_scoreboard_pkg::current_state
add wave -noupdate -expand -group ref_Signals -radix unsigned /count_scoreboard_pkg::timer_counter
add wave -noupdate -expand -group ref_Signals -radix unsigned /count_scoreboard_pkg::flag_counter
add wave -noupdate -expand -group ref_Signals -radix unsigned /count_scoreboard_pkg::count_value_ref
add wave -noupdate -expand -group ref_Signals -radix unsigned /count_scoreboard_pkg::busy_ref
TreeUpdate [SetDefaultTree]
WaveRestoreCursors {{Cursor 1} {70010037 ps} 0}
quietly wave cursor active 1
configure wave -namecolwidth 329
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
WaveRestoreZoom {69987968 ps} {70034804 ps}
