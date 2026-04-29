`timescale 1ns / 1ps

module TX_TOP #(
    parameter WL    = 16,         // System Word Length (Unscaled Q8.8)
    parameter N     = 4096,       // FFT Points
    parameter DDS_W = 8           // DDS resolution
)(
    input  wire                 clk,
    input  wire                 rst_n,
    
    // DDS Controls
    input  wire                 dds_enable,
    input  wire [31:0]          FTW_start,
    input  wire [12:0]          cycles,
    input  wire [31:0]          FTW_step,
    
    // OFDM Data Input (16-bit)
    input  wire signed [WL-1:0] ofdm_in_re,
    input  wire signed [WL-1:0] ofdm_in_im,
    output wire                 ofdm_rd_en,
    
    // Final TX Output (Time Domain)
    output wire                 tx_valid,
    output wire signed [WL-1:0] tx_out_re,
    output wire signed [WL-1:0] tx_out_im
);

    // Interconnect Wires
    wire                 dds_valid;
    wire [DDS_W-1:0]     dds_amplitude;
    wire                 fft_valid;
    wire signed [WL-1:0] fft_re, fft_im;
    wire                 bit_rev_valid;
    wire signed [WL-1:0] bit_rev_re, bit_rev_im;
    wire                 mux_valid;
    wire signed [WL-1:0] mux_re, mux_im;

    // =======================================================
    // 1. DDS Chirp Source
    // =======================================================
    dds_top #(
        .MEMORY_WIDTH(DDS_W)
    ) u_dds (
        .clk(clk), .rst_n(rst_n), 
        .enable(dds_enable),
        .FTW_start(FTW_start), .cycles(cycles), .FTW_step(FTW_step),
        .valid_out(dds_valid),
        .final_amplitude(dds_amplitude)
    );


    // Synchronized FFT Inputs
    wire                 fft_valid_in = fft_active;
    // Pass data when DDS is valid. Pass 16'd0 when DDS finishes early.
    wire signed [WL-1:0] fft_in_re    = (dds_valid) ? { {8{dds_amplitude[7]}}, dds_amplitude } : 16'd0;
    wire signed [WL-1:0] fft_in_im    = 16'd0;

    // =======================================================
    // 2. TX FFT (Analyzes Zero-Padded Chirp)
    // =======================================================
    fft_4096_top #(.WL(WL)) u_fft_tx (
        .clk(clk), .rst_n(rst_n),
        .valid_in(dds_valid),   // Driven by wrapper
        .in_real(fft_in_re),       // Driven by wrapper
        .in_imag(fft_in_im),       // Driven by wrapper
        .valid_out(fft_valid),
        .out_real(fft_re), .out_imag(fft_im)
    );

    // =======================================================
    // 3. Bit Reversal (Restore Frequency Order)
    // =======================================================
    bit_reversal_pingpong #(.WL(WL), .N(N)) u_bit_rev (
        .clk(clk), .rst_n(rst_n),
        .valid_in(fft_valid),
        .in_real(fft_re), .in_imag(fft_im),
        .valid_out(bit_rev_valid),
        .out_real(bit_rev_re), .out_imag(bit_rev_im)
    );

    // =======================================================
    // 4. ISAC MUX (Combine Radar & OFDM)
    // =======================================================
    isac_mux #(.WL(WL), .N(N)) u_mux (
        .clk(clk), .rst_n(rst_n),
        .radar_valid(bit_rev_valid),
        .radar_in_re(bit_rev_re>>>9), .radar_in_im(bit_rev_im>>>9),
        .ofdm_in_re(ofdm_in_re), .ofdm_in_im(ofdm_in_im),
        .ofdm_rd_en(ofdm_rd_en),
        .mux_valid(mux_valid),
        .mux_out_re(mux_re), .mux_out_im(mux_im)
    );

    // =======================================================
    // 5. TX IFFT (Transmit Signal Generation)
    // =======================================================
    ifft_4096_top #(.WL(WL)) u_ifft_tx (
        .clk(clk), .rst_n(rst_n),
        .valid_in(mux_valid),
        .in_real(mux_re), .in_imag(mux_im),
        .valid_out(tx_valid),
        .out_real(tx_out_re), .out_imag(tx_out_im)
    );

endmodule
