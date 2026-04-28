`timescale 1ns / 1ps

module isac_mux #(
    parameter WL = 16,
    parameter N = 4096
)(
    input  wire                 clk,
    input  wire                 rst_n,
    
    input  wire                 radar_valid,
    input  wire signed [WL-1:0] radar_in_re,
    input  wire signed [WL-1:0] radar_in_im,
    
    input  wire signed [WL-1:0] ofdm_in_re,
    input  wire signed [WL-1:0] ofdm_in_im,
    output reg                  ofdm_rd_en,
    
    output reg                  mux_valid,
    output reg  signed [WL-1:0] mux_out_re,
    output reg  signed [WL-1:0] mux_out_im
);

    // Buffer to hold radar data while OFDM is being transmitted
    reg signed [WL-1:0] radar_ram_re [0:2047];
    reg signed [WL-1:0] radar_ram_im [0:2047];
    reg [11:0] wr_ptr;
    
    always @(posedge clk) begin
        if (radar_valid) begin
            radar_ram_re[wr_ptr] <= radar_in_re;
            radar_ram_im[wr_ptr] <= radar_in_im;
            wr_ptr <= wr_ptr + 1;
        end else begin
            wr_ptr <= 0;
        end
    end

    // State Machine triggered by radar_valid rising edge
    reg [1:0] state;
    reg [11:0] rd_cnt;
    reg radar_valid_d;
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= 0;
            rd_cnt <= 0;
            ofdm_rd_en <= 0;
            mux_valid <= 0;
            radar_valid_d <= 0;
        end else begin
            radar_valid_d <= radar_valid;
            
            case (state)
                0: begin
                    mux_valid <= 0;
                    ofdm_rd_en <= 0;
                    if (radar_valid && !radar_valid_d) begin
                        state <= 1;
                        rd_cnt <= 0;
                        ofdm_rd_en <= 1; 
                    end
                end
                1: begin // Output 2048 OFDM bins
                    mux_valid <= 1;
                    mux_out_re <= ofdm_in_re;
                    mux_out_im <= ofdm_in_im;
                    rd_cnt <= rd_cnt + 1;
                    if (rd_cnt == 2046) ofdm_rd_en <= 0;
                    if (rd_cnt == 2047) begin
                        state <= 2;
                        rd_cnt <= 0;
                    end
                end
                2: begin // Output 381 Zeros (Padding)
                    mux_valid <= 1;
                    mux_out_re <= 0;
                    mux_out_im <= 0;
                    rd_cnt <= rd_cnt + 1;
                    if (rd_cnt == 380) begin
                        state <= 3;
                        rd_cnt <= 0;
                    end
                end
                3: begin // Output 1667 Radar bins
                    mux_valid <= 1;
                    mux_out_re <= radar_ram_re[rd_cnt];
                    mux_out_im <= radar_ram_im[rd_cnt];
                    rd_cnt <= rd_cnt + 1;
                    if (rd_cnt == 1667) begin
                        state <= 0;
                        mux_valid <= 0;
                    end
                end
            endcase
        end
    end
endmodule