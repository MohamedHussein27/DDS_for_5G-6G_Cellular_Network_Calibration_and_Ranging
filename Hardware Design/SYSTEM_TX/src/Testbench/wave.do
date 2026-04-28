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
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/clk
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/rst_n
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/valid_in
add wave -noupdate -group FFT -radix unsigned /TX_TOP_tb/uut/u_fft_tx/in_real
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/in_imag
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/valid_out
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/out_real
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/out_imag
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/global_addr
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/global_sel
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/pipeline_en
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/stage_re
add wave -noupdate -group FFT /TX_TOP_tb/uut/u_fft_tx/stage_im
add wave -noupdate /TX_TOP_tb/uut/u_fft_tx/out_real
add wave -noupdate -radix unsigned /TX_TOP_tb/uut/u_dds/u_phase_acc/cycles_counter
add wave -noupdate /TX_TOP_tb/uut/u_fft_tx/stg1/bf_inst/a_real
add wave -noupdate -radix unsigned /TX_TOP_tb/uut/u_fft_tx/ctrl/count
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/clk
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/rst_n
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/valid_in
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/in_real
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/in_imag
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/valid_out
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/out_real
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/out_imag
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/ram_ping_re
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/ram_ping_im
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/ram_pong_re
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/ram_pong_im
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/wr_ptr
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/rd_ptr
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/bank_sel
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/rd_bank
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/reading
add wave -noupdate -group Bit_Rev /TX_TOP_tb/uut/u_bit_rev/rd_addr
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/clk
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/rst_n
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/radar_valid
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/radar_in_re
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/radar_in_im
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/ofdm_in_re
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/ofdm_in_im
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/ofdm_rd_en
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/mux_valid
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/mux_out_re
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/mux_out_im
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/radar_ram_re
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/radar_ram_im
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/wr_ptr
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/state
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/rd_cnt
add wave -noupdate -group MUX /TX_TOP_tb/uut/u_mux/radar_valid_d
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/clk
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/rst_n
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/valid_in
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/in_real
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/in_imag
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/valid_out
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/out_real
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/out_imag
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/global_addr
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/global_sel
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/pipeline_en
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/stage_re
add wave -noupdate -expand -group IFFT /TX_TOP_tb/uut/u_ifft_tx/stage_im
TreeUpdate [SetDefaultTree]
WaveRestoreCursors {{Cursor 1} {24607985 ps} 0}
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
WaveRestoreZoom {24602631 ps} {24634631 ps}
