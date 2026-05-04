`timescale 1ns / 1ps

module TX_TOP #(
    parameter WL    = 16,
    parameter N     = 4096,
    parameter DDS_W = 8
)(
    input  wire                 clk,
    input  wire                 rst_n,
    input  wire [1:0]           dds_count,
    
    // Final TX Output (Natural-Order Time Domain)
    output wire                 tx_valid,
    output wire signed [WL-1:0] tx_out_re,
    output wire signed [WL-1:0] tx_out_im,

    // Bit-Reversal Outputs (Reference to RX)
    output wire                 bit_rev_valid_out,
    output wire signed [WL-1:0] bit_rev_out_re,
    output wire signed [WL-1:0] bit_rev_out_im
);

    reg [31:0]          FTW_start;
    reg [12:0]          cycles;
    reg [31:0]          FTW_step;
    
    // Changed to reg to fix multiple-driver error during synthesis
    reg                 dds_enable; 
    
    // Added 1-cycle delay to align with ROM latency
    reg [1:0]           dds_count_d; 
    wire [31:0]         config_data;

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

    wire                 ofdm_rd_en;
    wire signed [WL-1:0] ofdm_in_re, ofdm_in_im;
    reg  [10:0]          ofdm_ptr;

    // Address alignment pipeline for ROM
    always @(posedge clk or negedge rst_n) begin
        if(!rst_n) dds_count_d <= 2'b00;
        else       dds_count_d <= dds_count;
    end

    // Configuration FSM using aligned count
    always @(posedge clk or negedge rst_n) begin
        if(!rst_n) begin
            FTW_start  <= 32'd0;
            cycles     <= 13'd0;
            FTW_step   <= 32'd0;
            dds_enable <= 1'b0;
        end else begin
            case (dds_count_d) // Check against the delayed count
                2'b00: FTW_start  <= config_data; 
                2'b01: cycles     <= config_data[12:0];
                2'b10: FTW_step   <= config_data; 
                2'b11: dds_enable <= 1'b1; 
                default: begin
                    FTW_start <= 32'd0;
                    cycles    <= 13'd0;
                    FTW_step  <= 32'd0;
                end
            endcase
        end
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            ofdm_ptr <= 0;
        end else if (ofdm_rd_en) begin
            ofdm_ptr <= ofdm_ptr + 1;
        end
    end

    Config_reg u_config_reg (
        .addra(dds_count), // Address goes in immediately
        .clka(clk),
        .douta(config_data) // Data comes out one cycle later
    );

    OFDM_real u_ofdm_rom_real (
        .clka(clk),
        .addra(ofdm_ptr),
        .douta(ofdm_in_re)
    );

    OFDM_imag u_ofdm_rom_imag (
        .clka(clk),
        .addra(ofdm_ptr),
        .douta(ofdm_in_im)
    );

    // =======================================================
    // 1. DDS Chirp Source
    // =======================================================
    dds_top #(
        .TUNING_WORD_WIDTH(32), 
        .CYCLES_WIDTH(13),
        .ADDRESS_WIDTH(16), 
        .MEMORY_WIDTH(DDS_W)
    ) u_dds (
        .clk(clk), 
        .rst_n(rst_n), 
        .enable(dds_enable),
        .FTW_start(FTW_start), 
        .cycles(cycles), 
        .FTW_step(FTW_step),
        .valid_out(dds_valid),
        .final_amplitude(dds_amplitude)
    );

    // Sign extend the DDS amplitude to 16 bits
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
            mux_cnt <= 12'd0; 
        end
    end

    assign bit_rev_valid_out = mux_valid && (mux_cnt >= 2048);
    assign bit_rev_out_re    = mux_re;
    assign bit_rev_out_im    = mux_im;

endmodule