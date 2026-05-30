`timescale 1ns / 1ps

module SYSTEM_TOP #(
    parameter WL    = 16,
    parameter N     = 4096,
    parameter DDS_W = 8
)(
    input  wire                 clk,
    input  wire                 rst_n,
    
    // Master System Controls
    input  wire                 dds_enable,
    input  wire [31:0]          FTW_start,
    input  wire [12:0]          cycles,
    input  wire [31:0]          FTW_step,
    
    // Note: OFDM Data ports removed because ROM is now inside TX_TOP
    
    // Receiver Output: Communication Data
    output wire                 ofdm_valid_out,
    output wire signed [WL-1:0] ofdm_out_re,
    output wire signed [WL-1:0] ofdm_out_im,
    
    // Receiver Output: Radar Range Profile
    output wire                 radar_valid_out,
    output wire signed [WL-1:0] radar_out_re,
    output wire signed [WL-1:0] radar_out_im
);

    // =======================================================
    // 1. INTERCONNECT WIRES
    // =======================================================
    
    // Main Over-The-Air RF Interface (TX Output -> RX Input)
    wire                 rf_valid;
    wire signed [WL-1:0] rf_re;
    wire signed [WL-1:0] rf_im;
    
    // Internal DSP Reference Interface (TX DSP -> RX DSP)
    wire                 ref_valid;
    wire signed [WL-1:0] ref_re;
    wire signed [WL-1:0] ref_im;

    // RX Scaled Input
    wire signed [WL-1:0] rx_in_re;
    wire signed [WL-1:0] rx_in_im;

    // =======================================================
    // 2. TRANSMITTER (TX_TOP)
    // =======================================================
    TX_TOP #(
        .WL(WL), .N(N), .DDS_W(DDS_W)
    ) u_tx (
        .clk(clk), .rst_n(rst_n),
        .dds_enable(dds_enable),
        .FTW_start(FTW_start), .cycles(cycles), .FTW_step(FTW_step),
        
        // Output to RX internal reference
        .bit_rev_valid_out(ref_valid),
        .bit_rev_out_re(ref_re),
        .bit_rev_out_im(ref_im),
        
        // Output to Ideal Channel (Wire)
        .tx_valid(rf_valid),
        .tx_out_re(rf_re), .tx_out_im(rf_im)
    );

    // =======================================================
    // 3. THE CHANNEL LOOPBACK & SCALING (Q8.8 -> Q2.14)
    // =======================================================
    // Applies the division by 4096 and format conversion using a 6-bit 
    // arithmetic right shift. We add 32 (2^5) to round to the nearest integer 
    // rather than truncating, which preserves SQNR.
    
    assign rx_in_re = ($signed(rf_re) + 16'sd32) >>> 6;
    assign rx_in_im = ($signed(rf_im) + 16'sd32) >>> 6;

    // =======================================================
    // 4. RECEIVER (RX_TOP)
    // =======================================================
    RX_TOP #(
        .WL(WL),
        .N(N)
    ) u_rx (
        .clk(clk),
        .rst_n(rst_n),
        
        // Channel Input (Using the Scaled Wires!)
        .rx_valid_in(rf_valid),
        .rx_in_re(rx_in_re),
        .rx_in_im(rx_in_im),
        
        // Reference Input
        .ref_wr_en(ref_valid),
        .ref_wr_re(ref_re),
        .ref_wr_im(ref_im),
        
        // Outputs
        .ofdm_valid_out(ofdm_valid_out),
        .ofdm_out_re(ofdm_out_re),
        .ofdm_out_im(ofdm_out_im),
        
        .radar_valid_out(radar_valid_out),
        .radar_out_re(radar_out_re),
        .radar_out_im(radar_out_im)
    );

endmodule