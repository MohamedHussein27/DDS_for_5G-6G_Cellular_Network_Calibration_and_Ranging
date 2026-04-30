`timescale 1ns / 1ps

module isac_demux #(parameter WL = 16) (
    input  wire                 clk,
    input  wire                 rst_n,
    input  wire                 valid_in,
    input  wire signed [WL-1:0] in_re,
    input  wire signed [WL-1:0] in_im,
    
    output reg                  ofdm_valid,
    output reg  signed [WL-1:0] ofdm_out_re,
    output reg  signed [WL-1:0] ofdm_out_im,
    
    output reg                  radar_valid,
    output reg  signed [WL-1:0] radar_out_re,
    output reg  signed [WL-1:0] radar_out_im
);

    reg [11:0] count;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            count        <= 0;
            ofdm_valid   <= 0;
            radar_valid  <= 0;
            ofdm_out_re  <= 0;
            ofdm_out_im  <= 0;
            radar_out_re <= 0;
            radar_out_im <= 0;
        end else if (valid_in) begin
            count <= count + 1;
            if (count < 2048) begin
                ofdm_valid   <= 1'b1;
                ofdm_out_re  <= in_re;
                ofdm_out_im  <= in_im;
                radar_valid  <= 1'b0;
            end else begin
                radar_valid  <= 1'b1;
                radar_out_re <= in_re;
                radar_out_im <= in_im;
                ofdm_valid   <= 1'b0;
            end
        end else begin
            ofdm_valid  <= 0;
            radar_valid <= 0;
            count <= 0;
        end
    end
endmodule