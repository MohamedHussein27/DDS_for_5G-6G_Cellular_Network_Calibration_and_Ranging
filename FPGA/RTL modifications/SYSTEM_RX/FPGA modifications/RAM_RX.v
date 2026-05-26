`timescale 1ns / 1ps

module ref_ram_2048 #(
    parameter WL = 16
)(
    input  wire                 clk,
    input  wire                 rst_n,      // <-- ADDED Reset
    
    // --- PORT A: Write Interface (From TX) ---
    input  wire                 wr_en,
    input  wire signed [WL-1:0] wr_re,
    input  wire signed [WL-1:0] wr_im,
    
    // --- PORT B: Read Interface (To RX Multiplier) ---
    input  wire                 rd_en,      // <-- ADDED Read Enable
    output reg  signed [WL-1:0] rd_re,
    output reg  signed [WL-1:0] rd_im
);

    // Memory arrays for Real and Imaginary components
    reg signed [WL-1:0] ram_re [0:2047];
    reg signed [WL-1:0] ram_im [0:2047];

    // ==========================================
    // PORT A: Smart Write Logic (First 2048 Only)
    // ==========================================
    reg [11:0] wr_counter; 
    
    always @(posedge clk) begin
        if (!rst_n) begin
            wr_counter <= 12'd0;
        end else begin
            if (wr_en) begin
                // Only write if we haven't reached 2048 yet
                if (wr_counter < 12'd2048) begin
                    ram_re[wr_counter[10:0]] <= wr_re;
                    ram_im[wr_counter[10:0]] <= wr_im;
                    wr_counter <= wr_counter + 1;
                end
            end else begin
                // Reset counter between frames
                wr_counter <= 12'd0;
            end
        end
    end

    // ==========================================
    // PORT B: Smart Read Logic (Linear Read)
    // ==========================================
    reg [10:0] rd_counter;

    always @(posedge clk) begin
        if (!rst_n) begin
            rd_counter <= 11'd0;
            rd_re      <= 0;
            rd_im      <= 0;
        end else begin
            if (rd_en) begin
                rd_re      <= ram_re[rd_counter];
                rd_im      <= ram_im[rd_counter];
                rd_counter <= rd_counter + 1;
            end else begin
                // Reset counter between frames
                rd_counter <= 11'd0;
            end
        end
    end

endmodule