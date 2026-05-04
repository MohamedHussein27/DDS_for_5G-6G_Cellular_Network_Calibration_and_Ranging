`timescale 1ns / 1ps

module TX_wrapper #(
    parameter WL    = 16,
    parameter N     = 4096,
    parameter DDS_W = 8
)(
    input  wire                 clk,
    input  wire                 rst_n,
    
    // DDS Controls
    input  wire                 dds_enable,

    // Final Comparator Status Flags (Map these to physical LEDs!)
    output wire                 tx_valid,
    output wire                 real_match_flag,
    output wire                 imag_match_flag
);



 wire signed [WL-1:0] tx_out_re;
 wire signed [WL-1:0] tx_out_im;
    wire signed [WL-1:0] tx_ref_re_out;
    wire signed [WL-1:0] tx_ref_im_out;
    wire enable_out;
    wire [1:0] dds_count;

    wire [11:0] write_addr; 
    wire signed [WL-1:0] tx_ref_re;
    wire signed [WL-1:0] tx_ref_im;
(* dont_touch = "true" *)    wire bit_rev_valid_out;
(* dont_touch = "true" *)    wire signed [WL-1:0] bit_rev_out_re, bit_rev_out_im;

    // Map internal wires to output ports
    assign tx_ref_re_out = tx_ref_re;
    assign tx_ref_im_out = tx_ref_im;

    // ==============================================================
    // CRITICAL FIX: Pipeline registers to align TX data with 
    // the 1-clock-cycle read latency of the BRAM reference ROMs
    // ==============================================================
//    reg signed [WL-1:0] tx_out_re_d;
//    reg signed [WL-1:0] tx_out_im_d;
//    reg                 tx_valid_d;
    reg [11:0]          write_addr_d;

     always @(posedge clk or negedge rst_n) begin
         if (!rst_n) begin
    //         tx_out_re_d  <= 0;
    //         tx_out_im_d  <= 0;
    //         tx_valid_d   <= 0;
             write_addr_d <= 0;
         end else begin
    //         tx_out_re_d  <= tx_out_re;
    //         tx_out_im_d  <= tx_out_im;
    //         tx_valid_d   <= tx_valid;
             write_addr_d <= write_addr; // So comparator knows exact symbol
    //     end
     end
     end
    // --------------------------------------------------------------
    // Instantiations
    // --------------------------------------------------------------
    dds_counter #(
        .ADDR_WIDTH(2)
    ) u_dds_counter (
        .enable(dds_enable),
        .clk(clk),
        .rst_n(rst_n),
        .count(dds_count)
//        .enable_out(enable_out)
    );

    TX_TOP #(
        .WL(WL),
        .N(N),
        .DDS_W(DDS_W)
    ) u_tx_top (
        .clk(clk),
        .rst_n(rst_n),
        .dds_count(dds_count),
        .tx_valid(tx_valid),
        .tx_out_re(tx_out_re),
        .tx_out_im(tx_out_im),
        .bit_rev_valid_out(bit_rev_valid_out),
        .bit_rev_out_re(bit_rev_out_re),
        .bit_rev_out_im(bit_rev_out_im)
    );

    bram_write_counter #(
        .ADDR_WIDTH(12)
    ) u_bram_write_counter (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(tx_valid),
        .write_addr(write_addr)
    );

    TX_ref_real u_ref_real_mem(
        .clka(clk),
        .addra(write_addr),
        .douta(tx_ref_re)
    );

    TX_ref_imag u_ref_imag_mem(
        .clka(clk),
        .addra(write_addr),
        .douta(tx_ref_im)
    );

    // Feed the DELAYED signals to perfectly match the ROM outputs
    real_comparator #(
        .INPUT_WIDTH(WL), 
        .NUM_SYMBOLS(N), 
        .SYMBOL(12) 
    ) u_real_comp (
        .clk(clk) , 
        .rst_n(rst_n),
        .valid_in(tx_valid),
        .fft_real_in(tx_out_re),
        .fft_real_ref_in(tx_ref_re),
        .symbol_index(write_addr_d),
        .comparison_result(real_match_flag)
    );

    imaginary_comparator #(
        .INPUT_WIDTH(WL), 
        .NUM_SYMBOLS(N), 
        .SYMBOL(12)
    ) u_imag_comp (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(tx_valid),
        .fft_imag_in(tx_out_im),
        .fft_imag_ref_in(tx_ref_im),
        .symbol_index(write_addr_d),
        .comparison_result(imag_match_flag)
    );

endmodule