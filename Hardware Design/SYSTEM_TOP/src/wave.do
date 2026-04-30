onerror {resume}
quietly WaveActivateNextPane {} 0
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/clk
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/rst_n
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/dds_enable
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/FTW_start
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/cycles
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/FTW_step
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/ofdm_in_re
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/ofdm_in_im
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/ofdm_rd_en
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/tx_valid
add wave -noupdate -group TX -radix decimal /SYSTEM_TOP_tb/uut/u_tx/tx_out_re
add wave -noupdate -group TX -radix decimal /SYSTEM_TOP_tb/uut/u_tx/tx_out_im
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/bit_rev_valid_out
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/bit_rev_out_re
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/bit_rev_out_im
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/dds_valid
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/dds_amplitude
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/fft_valid
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/fft_re
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/fft_im
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/bit_rev_valid
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/bit_rev_re
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/bit_rev_im
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/mux_valid
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/mux_re
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/mux_im
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/fft_in_re
add wave -noupdate -group TX /SYSTEM_TOP_tb/uut/u_tx/fft_in_im
add wave -noupdate -radix decimal /SYSTEM_TOP_tb/uut/u_tx/u_mux/mux_out_re
add wave -noupdate -radix decimal /SYSTEM_TOP_tb/uut/u_tx/u_mux/mux_out_im
add wave -noupdate /SYSTEM_TOP_tb/uut/u_tx/u_mux/mux_valid
add wave -noupdate -expand -group RX_RAM /SYSTEM_TOP_tb/uut/u_rx_top/u_ref_ram/clk
add wave -noupdate -expand -group RX_RAM /SYSTEM_TOP_tb/uut/u_rx_top/u_ref_ram/wr_en
add wave -noupdate -expand -group RX_RAM /SYSTEM_TOP_tb/uut/u_rx_top/u_ref_ram/wr_re
add wave -noupdate -expand -group RX_RAM /SYSTEM_TOP_tb/uut/u_rx_top/u_ref_ram/wr_im
add wave -noupdate -expand -group RX_RAM /SYSTEM_TOP_tb/uut/u_rx_top/u_ref_ram/rd_re
add wave -noupdate -expand -group RX_RAM /SYSTEM_TOP_tb/uut/u_rx_top/u_ref_ram/rd_im
add wave -noupdate -expand -group RX_RAM /SYSTEM_TOP_tb/uut/u_rx_top/u_ref_ram/ram_re
add wave -noupdate -expand -group RX_RAM /SYSTEM_TOP_tb/uut/u_rx_top/u_ref_ram/ram_im
add wave -noupdate -group RX_MULT -radix hexadecimal /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/re1
add wave -noupdate -group RX_MULT /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/im1
add wave -noupdate -group RX_MULT -radix hexadecimal /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/re2
add wave -noupdate -group RX_MULT /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/im2
add wave -noupdate -group RX_MULT -radix hexadecimal /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/re_out
add wave -noupdate -group RX_MULT /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/im_out
add wave -noupdate -group RX_MULT /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/mul1
add wave -noupdate -group RX_MULT /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/mul2
add wave -noupdate -group RX_MULT /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/mul3
add wave -noupdate -group RX_MULT /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/mul4
add wave -noupdate -group RX_MULT /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/sub_re
add wave -noupdate -group RX_MULT /SYSTEM_TOP_tb/uut/u_rx_top/u_mult/add_im
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/clk
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/rst_n
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/valid_in
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/in_re
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/in_im
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/ofdm_valid
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/ofdm_out_re
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/ofdm_out_im
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/radar_valid
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/radar_out_re
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/radar_out_im
add wave -noupdate -group RX_demux /SYSTEM_TOP_tb/uut/u_rx_top/u_demux/count
add wave -noupdate /SYSTEM_TOP_tb/uut/u_rx_top/mult_valid
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/clk
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/rst_n
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/valid_in
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/in_real
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/in_imag
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/valid_out
add wave -noupdate -expand -group RX_IFFT -radix decimal /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/out_real
add wave -noupdate -expand -group RX_IFFT -radix decimal /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/out_imag
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/global_addr
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/global_sel
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/pipeline_en
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/stage_re
add wave -noupdate -expand -group RX_IFFT /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_ifft/stage_im
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/clk
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/rst_n
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/valid_in
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/in_real
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/in_imag
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/valid_out
add wave -noupdate -expand -group RX_REV -radix decimal /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/out_real
add wave -noupdate -expand -group RX_REV -radix decimal /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/out_imag
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/wr_ptr
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/rd_ptr
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/bank_sel
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/rd_bank
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/reading
add wave -noupdate -expand -group RX_REV /SYSTEM_TOP_tb/uut/u_rx_top/u_rx_final_rev/rd_addr
TreeUpdate [SetDefaultTree]
WaveRestoreCursors {{Cursor 1} {45107837 ps} 0}
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
WaveRestoreZoom {45097347 ps} {45129347 ps}
