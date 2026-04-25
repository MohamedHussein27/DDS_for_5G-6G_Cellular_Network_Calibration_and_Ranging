onerror {resume}
quietly WaveActivateNextPane {} 0
add wave -noupdate -group TX /TX_TOP_tb/clk
add wave -noupdate -group TX /TX_TOP_tb/rst_n
add wave -noupdate -group TX /TX_TOP_tb/dds_enable
add wave -noupdate -group TX /TX_TOP_tb/FTW_start
add wave -noupdate -group TX /TX_TOP_tb/cycles
add wave -noupdate -group TX /TX_TOP_tb/FTW_step
add wave -noupdate -group TX /TX_TOP_tb/ofdm_in_re
add wave -noupdate -group TX /TX_TOP_tb/ofdm_in_im
add wave -noupdate -group TX /TX_TOP_tb/ofdm_rd_en
add wave -noupdate -group TX /TX_TOP_tb/tx_valid
add wave -noupdate -group TX /TX_TOP_tb/tx_out_re
add wave -noupdate -group TX /TX_TOP_tb/tx_out_im
add wave -noupdate -group TX /TX_TOP_tb/ofdm_rom_re
add wave -noupdate -group TX /TX_TOP_tb/ofdm_rom_im
add wave -noupdate -group TX /TX_TOP_tb/ofdm_ptr
add wave -noupdate -group TX /TX_TOP_tb/file_re
add wave -noupdate -group TX /TX_TOP_tb/file_im
add wave -noupdate -group TX /TX_TOP_tb/sample_count
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/clk
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/rst_n
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/enable
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/FTW_start
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/cycles
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/FTW_step
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/valid_out
add wave -noupdate -expand -group DDS -radix unsigned /TX_TOP_tb/uut/u_dds/final_amplitude
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/phase_address
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/mapped_address
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/neg_flag
add wave -noupdate -expand -group DDS -radix decimal /TX_TOP_tb/uut/u_dds/amplitude_out
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/valid_phase
add wave -noupdate -expand -group DDS /TX_TOP_tb/uut/u_dds/valid_out_reg
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/clk
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/rst_n
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/valid_in
add wave -noupdate -expand -group FFT -radix unsigned /TX_TOP_tb/uut/u_fft_tx/in_real
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/in_imag
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/valid_out
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/out_real
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/out_imag
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/global_addr
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/global_sel
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/pipeline_en
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/stage_re
add wave -noupdate -expand -group FFT /TX_TOP_tb/uut/u_fft_tx/stage_im
add wave -noupdate /TX_TOP_tb/uut/u_fft_tx/out_real
add wave -noupdate -radix unsigned /TX_TOP_tb/uut/u_dds/u_phase_acc/cycles_counter
add wave -noupdate /TX_TOP_tb/uut/u_fft_tx/stg1/bf_inst/a_real
add wave -noupdate /TX_TOP_tb/uut/u_fft_tx/stg1/bf_inst/enable
add wave -noupdate -radix unsigned /TX_TOP_tb/uut/u_fft_tx/ctrl/count
add wave -noupdate /TX_TOP_tb/uut/u_fft_tx/ctrl/count_valid
TreeUpdate [SetDefaultTree]
WaveRestoreCursors {{Cursor 1} {8218483 ps} 0}
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
configure wave -timelineunits ps
update
WaveRestoreZoom {8209717 ps} {8241717 ps}
