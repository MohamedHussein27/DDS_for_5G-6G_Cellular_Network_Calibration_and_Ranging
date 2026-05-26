// `timescale 1ns / 1ps

// module RX_TOP #(
//     parameter WL = 16,
//     parameter N  = 4096 // Keep this matching your system size
// )(
//     input  wire                 clk,
//     input  wire                 rst_n,
    
//     // 1. Input from Channel (Time Domain from TX)
//     input  wire                 rx_valid_in,
//     input  wire signed [WL-1:0] rx_in_re,
//     input  wire signed [WL-1:0] rx_in_im,
    
//     // 2. Dynamic Reference RAM Interface (From TX)
//     input  wire                 ref_wr_en,
//     input  wire signed [WL-1:0] ref_wr_re,
//     input  wire signed [WL-1:0] ref_wr_im,
    
//     // 3. OFDM Output (Communication)
//     output wire                 ofdm_valid_out,
//     output wire signed [WL-1:0] ofdm_out_re,
//     output wire signed [WL-1:0] ofdm_out_im,
    
//     // 4. Radar Output (Range Profile)
//     output wire                 radar_valid_out,
//     output wire signed [WL-1:0] radar_out_re,
//     output wire signed [WL-1:0] radar_out_im
// );

//     // =======================================================
//     // 1. RX FFT (Forward Transform)
//     // =======================================================
//     wire                 fft_valid;
//     wire signed [WL-1:0] fft_re, fft_im;
    
//     fft_4096_top #(.WL(WL)) u_rx_fft (
//         .clk(clk), 
//         .rst_n(rst_n),
//         .valid_in(rx_valid_in),
//         .in_real(rx_in_re), 
//         .in_imag(rx_in_im),
//         .valid_out(fft_valid),
//         .out_real(fft_re), 
//         .out_imag(fft_im)
//     );

//     // =======================================================
//     // 2. ISAC Demultiplexer
//     // =======================================================
//     wire                 radar_split_valid;
//     wire signed [WL-1:0] radar_split_re, radar_split_im;

//     isac_demux #(.WL(WL)) u_demux (
//         .clk(clk), 
//         .rst_n(rst_n),
//         .valid_in(fft_valid),
//         .in_re(fft_re), 
//         .in_im(fft_im),
        
//         .ofdm_valid(ofdm_valid_out),
//         .ofdm_out_re(ofdm_out_re), 
//         .ofdm_out_im(ofdm_out_im),
        
//         .radar_valid(radar_split_valid),
//         .radar_out_re(radar_split_re), 
//         .radar_out_im(radar_split_im)
//     );

//     // =======================================================
//     // 3. REFERENCE RAM
//     // =======================================================
//     wire signed [WL-1:0] ram_out_re, ram_out_im;

//     ref_ram_2048 #(.WL(WL)) u_ref_ram (
//         .clk(clk),
//         .rst_n(rst_n),
        
//         // Write Port (From TX)
//         .wr_en(ref_wr_en),
//         .wr_re(ref_wr_re),
//         .wr_im(ref_wr_im),
        
//         // Read Port (Triggered by radar split valid)
//         .rd_en(radar_split_valid),
//         .rd_re(ram_out_re),
//         .rd_im(ram_out_im)
//     );

//     // =======================================================
//     // 4. DELAY MATCHING & CONJUGATE MULTIPLIER
//     // =======================================================
//     // The RAM takes 1 cycle to output data, so we must delay the 
//     // radar data and the valid signal by 1 cycle to align with it.
//     reg signed [WL-1:0] radar_split_re_d, radar_split_im_d;
//     reg                 radar_split_valid_d;

//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             radar_split_valid_d <= 1'b0;
//             radar_split_re_d    <= 0;
//             radar_split_im_d    <= 0;
//         end else begin
//             radar_split_valid_d <= radar_split_valid;
//             radar_split_re_d    <= radar_split_re;
//             radar_split_im_d    <= radar_split_im;
//         end
//     end

//     // Perform the complex conjugate mathematically (invert imaginary)
//     wire signed [WL-1:0] ref_conj_im = -ram_out_im; 
//     wire signed [WL-1:0] mult_re, mult_im;
    
//     multiplier #(.WL(WL)) u_mult (
//         .re1(radar_split_re_d), .im1(radar_split_im_d),
//         .re2(ram_out_re>>>7),       .im2(ref_conj_im>>>7),
//         .re_out(mult_re),       .im_out(mult_im)
//     );

//     // The delayed valid signal acts as the valid output for the multiplier
//     wire mult_valid = radar_split_valid_d;

//     // =======================================================
//     // 5. RX IFFT (Inverse Transform for Range Profile)
//     // =======================================================
//     wire                 ifft_valid;
//     wire signed [WL-1:0] ifft_re, ifft_im;

//     ifft_2048_top #(.WL(WL)) u_rx_ifft (
//         .clk(clk), 
//         .rst_n(rst_n),
//         .valid_in(mult_valid),
//         .in_real(mult_re), 
//         .in_imag(mult_im),
//         .valid_out(ifft_valid),
//         .out_real(ifft_re), 
//         .out_imag(ifft_im)
//     );

//     // =======================================================
//     // 6. Final Bit Reversal (Natural Order Output)
//     // =======================================================
//     bit_reversal_pingpong #(.WL(WL), .N(N/2), .ADDR_W(11)) u_rx_final_rev (
//         .clk(clk), 
//         .rst_n(rst_n),
//         .valid_in(ifft_valid),
//         .in_real(ifft_re), 
//         .in_imag(ifft_im),
//         .valid_out(radar_valid_out),
//         .out_real(radar_out_re), 
//         .out_imag(radar_out_im)
//     );

// endmodule

`timescale 1ns / 1ps

module RX_TOP #(
    parameter WL = 16,
    parameter N  = 4096 // Keep this matching your system size
)(
    input  wire                 clk,
    input  wire                 rst_n,
    
    // 1. Input from Channel (Time Domain from TX)
    input  wire                 rx_valid_in,
    input  wire signed [WL-1:0] rx_in_re,
    input  wire signed [WL-1:0] rx_in_im,
    
    // 3. OFDM Output (Communication)
    output wire                 ofdm_valid_out,
    output wire signed [WL-1:0] ofdm_out_re,
    output wire signed [WL-1:0] ofdm_out_im,
    
    // 4. Radar Output (Range Profile)
    output wire                 radar_valid_out,
    output wire signed [WL-1:0] radar_out_re,
    output wire signed [WL-1:0] radar_out_im
);

    // =======================================================
    // 1. RX FFT (Forward Transform)
    // =======================================================
    wire                 fft_valid;
    wire signed [WL-1:0] fft_re, fft_im;
    
    fft_4096_top #(.WL(WL)) u_rx_fft (
        .clk(clk), 
        .rst_n(rst_n),
        .valid_in(rx_valid_in),
        .in_real(rx_in_re), 
        .in_imag(rx_in_im),
        .valid_out(fft_valid),
        .out_real(fft_re), 
        .out_imag(fft_im)
    );

    // =======================================================
    // NEW: 1.5 Mid-Path Bit Reversal (Unscramble Frequencies)
    // =======================================================
    wire                 mid_rev_valid;
    wire signed [WL-1:0] mid_rev_re, mid_rev_im;

    bit_reversal_pingpong #(
        .WL(WL), 
        .N(N), 
        .ADDR_W(12)
    ) u_rx_mid_rev (
        .clk(clk), 
        .rst_n(rst_n),
        .valid_in(fft_valid),
        .in_real(fft_re), 
        .in_imag(fft_im),
        .valid_out(mid_rev_valid),
        .out_real(mid_rev_re), 
        .out_imag(mid_rev_im)
    );

    // =======================================================
    // 2. ISAC Demultiplexer
    // =======================================================
    wire                 radar_split_valid;
    wire signed [WL-1:0] radar_split_re, radar_split_im;

    isac_demux #(.WL(WL)) u_demux (
        .clk(clk), 
        .rst_n(rst_n),
        .valid_in(mid_rev_valid),       // <--- Now fed by Bit Reversal
        .in_re(mid_rev_re),             // <--- Now fed by Bit Reversal
        .in_im(mid_rev_im),             // <--- Now fed by Bit Reversal
        
        .ofdm_valid(ofdm_valid_out),
        .ofdm_out_re(ofdm_out_re), 
        .ofdm_out_im(ofdm_out_im),
        
        .radar_valid(radar_split_valid),
        .radar_out_re(radar_split_re), 
        .radar_out_im(radar_split_im)
    );

    // =======================================================
    // 3. REFERENCE RAM
    // =======================================================
    wire signed [WL-1:0] ram_out_re, ram_out_im;

    // ref_ram_2048 #(.WL(WL)) u_ref_ram (
    //     .clk(clk),
    //     .rst_n(rst_n),
        
    //     // Write Port (From TX)
    //     .wr_en(ref_wr_en),
    //     .wr_re(ref_wr_re),
    //     .wr_im(ref_wr_im),
        
    //     // Read Port (Triggered by radar split valid)
    //     .rd_en(radar_split_valid),
    //     .rd_re(ram_out_re),
    //     .rd_im(ram_out_im)
    // );
    wire [10:0] ram_read_addr;
    rx_ram_counter #(.ADDR_WIDTH(11)) u_ref_ram_counter (
        .clk(clk),
        .rst_n(rst_n),
        .en(radar_split_valid),
        .count(ram_read_addr)
    );
    ref_ram_2048_re  u_ref_rom_re (
        .clka(clk),
        .addra(ram_read_addr),
        .douta(ram_out_re)
    );

    ref_ram_2048_im  u_ref_rom_im (
        .clka(clk),
        .addra(ram_read_addr),
        .douta(ram_out_im)
    );
    
    // =======================================================
    // 4. DELAY MATCHING & CONJUGATE MULTIPLIER
    // =======================================================
    // The RAM takes 1 cycle to output data, so we must delay the 
    // radar data and the valid signal by 1 cycle to align with it.
    reg signed [WL-1:0] radar_split_re_d, radar_split_im_d;
    reg                 radar_split_valid_d;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            radar_split_valid_d <= 1'b0;
            radar_split_re_d    <= 0;
            radar_split_im_d    <= 0;
        end else begin
            radar_split_valid_d <= radar_split_valid;
            radar_split_re_d    <= radar_split_re;
            radar_split_im_d    <= radar_split_im;
        end
    end

    // Perform the complex conjugate mathematically (invert imaginary)
    wire signed [WL-1:0] ref_conj_im = -ram_out_im; 
    wire signed [WL-1:0] mult_re, mult_im;
    
    multiplier_algorithm #(.WL(WL)) u_mult (
        .clk(clk),
        .rst_n(rst_n),
        .re1(radar_split_re_d), .im1(radar_split_im_d),
//        .re2(ram_out_re>>>7),       .im2(ref_conj_im>>>7),
        .re2(ram_out_re),       .im2(ref_conj_im),
        .re_out(mult_re),       .im_out(mult_im)
    );

    // The delayed valid signal acts as the valid output for the multiplier
    reg mult_valid;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            mult_valid <= 1'b0;
        end else begin
            mult_valid <= radar_split_valid_d;
        end
    end

    // =======================================================
    // 5. RX IFFT (Inverse Transform for Range Profile)
    // =======================================================
    wire                 ifft_valid;
    wire signed [WL-1:0] ifft_re, ifft_im;

    ifft_2048_top #(.WL(WL)) u_rx_ifft (
        .clk(clk), 
        .rst_n(rst_n),
        .valid_in(mult_valid),
        .in_real(mult_re>>>3), 
        .in_imag(mult_im>>>3),
//        .in_real(mult_re), 
//        .in_imag(mult_im),
        .valid_out(ifft_valid),
        .out_real(ifft_re), 
        .out_imag(ifft_im)
    );

    // =======================================================
    // 6. Final Bit Reversal (Natural Order Output)
    // =======================================================
    bit_reversal_pingpong #(.WL(WL), .N(N/2), .ADDR_W(11)) u_rx_final_rev (
        .clk(clk), 
        .rst_n(rst_n),
        .valid_in(ifft_valid),
        .in_real(ifft_re), 
        .in_imag(ifft_im),
        .valid_out(radar_valid_out),
        .out_real(radar_out_re), 
        .out_imag(radar_out_im)
    );

endmodule