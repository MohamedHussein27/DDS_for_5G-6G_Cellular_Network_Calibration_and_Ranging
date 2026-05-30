`timescale 1ns / 1ps

module TX_TOP #(
    parameter WL    = 16,
    parameter N     = 4096,
    parameter DDS_W = 8
)(
    input  wire                 clk,
    input  wire                 rst_n,
    
    // DDS Controls
    input  wire                 dds_enable,
    input  wire [31:0]          FTW_start,
    input  wire [12:0]          cycles,
    input  wire [31:0]          FTW_step,
    
    // Final TX Output (Natural-Order Time Domain)
    output wire                 tx_valid,
    output wire signed [WL-1:0] tx_out_re,
    output wire signed [WL-1:0] tx_out_im,

    // Bit-Reversal Outputs (Reference to RX)
    output wire                 bit_rev_valid_out,
    output wire signed [WL-1:0] bit_rev_out_re,
    output wire signed [WL-1:0] bit_rev_out_im
);

    // -------------------------------------------------------
    // Internal Wires
    // -------------------------------------------------------
    wire                 dds_valid;
    wire [DDS_W-1:0]     dds_amplitude;
    wire                 fft_valid;
    wire signed [WL-1:0] fft_re, fft_im;
    wire                 bit_rev_valid;
    wire signed [WL-1:0] bit_rev_re, bit_rev_im;
    wire                 mux_valid;
    wire signed [WL-1:0] mux_re, mux_im;
    wire                 ifft_valid;
    wire signed [WL-1:0] ifft_out_re, ifft_out_im;

    // =======================================================
    // NEW: Internal OFDM ROM & Pointer Logic
    // =======================================================
    wire                 ofdm_rd_en;
    wire signed [WL-1:0] ofdm_in_re, ofdm_in_im;
    reg [10:0]           ofdm_ptr;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            ofdm_ptr <= 0;
        end else if (ofdm_rd_en) begin
            ofdm_ptr <= ofdm_ptr + 1;
        end
    end

    ofdm_rom #(
        .DEPTH(2048),
        .WL(WL)
    ) u_ofdm_rom (
        .clk(clk),
        .addr(ofdm_ptr),
        .out_real(ofdm_in_re),
        .out_imag(ofdm_in_im)
    );

    // =======================================================
    // 1. DDS Chirp Source
    // =======================================================
    dds_top #(.MEMORY_WIDTH(DDS_W)) u_dds (
        .clk(clk), .rst_n(rst_n), 
        .enable(dds_enable),
        .FTW_start(FTW_start), .cycles(cycles), .FTW_step(FTW_step),
        .valid_out(dds_valid),
        .final_amplitude(dds_amplitude)
    );

    wire signed [WL-1:0] fft_in_re = (dds_valid) ? { {8{dds_amplitude[7]}}, dds_amplitude } : 16'd0;
    wire signed [WL-1:0] fft_in_im = 16'd0;

    // =======================================================
    // 2. TX FFT
    // =======================================================
    fft_4096_top #(.WL(WL)) u_fft_tx (
        .clk(clk), .rst_n(rst_n),
        .valid_in(dds_valid),
        .in_real(fft_in_re), .in_imag(fft_in_im),
        .valid_out(fft_valid),
        .out_real(fft_re), .out_imag(fft_im)
    );

    // =======================================================
    // 3. Frequency Domain Bit Reversal
    // =======================================================
    bit_reversal_pingpong #(.WL(WL), .N(N)) u_bit_rev (
        .clk(clk), .rst_n(rst_n),
        .valid_in(fft_valid),
        .in_real(fft_re), .in_imag(fft_im),
        .valid_out(bit_rev_valid),
        .out_real(bit_rev_re), .out_imag(bit_rev_im)
    );

    // =======================================================
    // 4. ISAC MUX
    // =======================================================
    isac_mux #(.WL(WL), .N(N)) u_mux (
        .clk(clk), .rst_n(rst_n),
        .radar_valid(bit_rev_valid),
        .radar_in_re(bit_rev_re >>> 7), .radar_in_im(bit_rev_im >>> 7),
        .ofdm_in_re(ofdm_in_re),        .ofdm_in_im(ofdm_in_im),
        .ofdm_rd_en(ofdm_rd_en),
        .mux_valid(mux_valid),
        .mux_out_re(mux_re),            .mux_out_im(mux_im)
    );

    // =======================================================
    // 5. TX IFFT
    // =======================================================
    ifft_4096_top #(.WL(WL)) u_ifft_tx (
        .clk(clk), .rst_n(rst_n),
        .valid_in(mux_valid),
        .in_real(mux_re), .in_imag(mux_im),
        .valid_out(ifft_valid),
        .out_real(ifft_out_re), .out_imag(ifft_out_im)
    );

    // =======================================================
    // 6. Time Domain Bit Reversal (Unscramble for DAC)
    // =======================================================
    bit_reversal_pingpong #(.WL(WL), .N(N)) u_tx_time_rev (
        .clk(clk), .rst_n(rst_n),
        .valid_in(ifft_valid),
        .in_real(ifft_out_re), .in_imag(ifft_out_im),
        .valid_out(tx_valid),
        .out_real(tx_out_re), .out_imag(tx_out_im)
    );

    // assign bit_rev_valid_out = bit_rev_valid;
    // assign bit_rev_out_re    = bit_rev_re;
    // assign bit_rev_out_im    = bit_rev_im;

    // =======================================================
    // Reference Output Qualification (Capture Last 2048)
    // =======================================================
    reg [11:0] mux_cnt;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            mux_cnt <= 12'd0;
        end else if (mux_valid) begin
            mux_cnt <= mux_cnt + 1;
        end else begin
            mux_cnt <= 12'd0; // Reset counter between 4096-bin frames
        end
    end

    // Only output 'valid' to the RX Reference RAM during the Radar half (2048-4095)
    assign bit_rev_valid_out = mux_valid && (mux_cnt >= 2048);
    assign bit_rev_out_re    = mux_re;
    assign bit_rev_out_im    = mux_im;

endmodule