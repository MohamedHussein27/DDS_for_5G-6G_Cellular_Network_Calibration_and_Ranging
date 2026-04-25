`timescale 1ns / 1ps

module isac_mux #(
    parameter WL = 16,
    parameter N  = 4096
)(
    input  wire                 clk,
    input  wire                 rst_n,
    
    // Radar Input (from Bit Reversal)
    input  wire                 radar_valid,
    input  wire signed [WL-1:0] radar_in_re,
    input  wire signed [WL-1:0] radar_in_im,
    
    // OFDM Input (from external source)
    input  wire signed [WL-1:0] ofdm_in_re,
    input  wire signed [WL-1:0] ofdm_in_im,
    output reg                  ofdm_rd_en,
    
    // Combined Output (to IFFT)
    output reg                  mux_valid,
    output reg  signed [WL-1:0] mux_out_re,
    output reg  signed [WL-1:0] mux_out_im
);

    // Internal Buffer for Radar Bins (First 2048 bins of FFT)
    reg signed [WL-1:0] radar_mem_re [0:2047];
    reg signed [WL-1:0] radar_mem_im [0:2047];
    reg [10:0]          wr_ptr;
    reg                 radar_ready;
    reg [11:0]          mux_cnt;

    // 1. Buffer Logic: Capture Radar samples
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            wr_ptr      <= 0;
            radar_ready <= 0;
        end else if (radar_valid && !radar_ready) begin
            radar_mem_re[wr_ptr] <= radar_in_re;
            radar_mem_im[wr_ptr] <= radar_in_im;
            wr_ptr <= wr_ptr + 1;
            if (wr_ptr == 2047) radar_ready <= 1'b1;
        end else if (mux_cnt == 4095 && mux_valid) begin
            radar_ready <= 0; // Reset for next frame
            wr_ptr      <= 0;
        end
    end

    // 2. Multiplexer Logic: Construct the 4096-bin frame
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            mux_cnt    <= 0;
            mux_valid  <= 0;
            ofdm_rd_en <= 0;
            mux_out_re <= 0;
            mux_out_im <= 0;
        end else if (radar_ready) begin
            mux_valid <= 1'b1;
            mux_cnt   <= mux_cnt + 1;
            
            if (mux_cnt < 2048) begin
                // Phase 1: OFDM Band (0-2047)
                // Shift right by 8 to match the radar guard-bit power level
                mux_out_re <= $signed(ofdm_in_re)>>>8;
                mux_out_im <= $signed(ofdm_in_im)>>>8;
                ofdm_rd_en <= 1'b1;
            end else begin
                // Phase 2: Radar Band (2048-4095)
                mux_out_re [9:0] <= radar_mem_re[mux_cnt[10:0]][15:6];
                mux_out_re [15:10]<= {6{radar_mem_re[15]}};
                mux_out_im <= radar_mem_im[mux_cnt[10:0]][15:6];
                mux_out_im <= {6{radar_mem_im[15]}};
                ofdm_rd_en <= 1'b0;
            end

            if (mux_cnt == 4095) mux_valid <= 1'b0;
        end else begin
            mux_valid  <= 1'b0;
            mux_cnt    <= 0;
            ofdm_rd_en <= 1'b0;
        end
    end

endmodule