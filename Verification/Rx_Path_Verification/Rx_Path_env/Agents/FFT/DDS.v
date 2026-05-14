`timescale 1ns / 1ps

module dds_top #(
    parameter TUNING_WORD_WIDTH = 32,
    parameter CYCLES_WIDTH      = 13,
    parameter ADDRESS_WIDTH     = 16, // Truncated phase width
    parameter MEMORY_WIDTH      = 8   // DAC resolution
)(
    input  wire                         clk,
    input  wire                         rst_n,
    input  wire                         enable,     // <--- NEW: Start the DDS chirp
    input  wire [TUNING_WORD_WIDTH-1:0] FTW_start,
    input  wire [CYCLES_WIDTH-1:0]      cycles,
    input  wire [TUNING_WORD_WIDTH-1:0] FTW_step,
    output wire                         valid_out,  // <--- NEW: Connects directly to FFT valid_in
    output wire signed [MEMORY_WIDTH-1:0]      final_amplitude
);

    // --- Internal Wires ---
    wire [ADDRESS_WIDTH-1:0] phase_address;
    wire [ADDRESS_WIDTH-1:0] mapped_address;
    wire                     neg_flag;
    wire [MEMORY_WIDTH-1:0]  amplitude_out;
    
    wire valid_phase;      // Raw valid signal from the phase accumulator
    reg  valid_out_reg;    // Pipeline register to match LUT latency

    // --- 1. Phase Accumulator ---
    PHASE_ACC #(
        .TUNING_WORD_WIDTH      (TUNING_WORD_WIDTH),
        .CYCLES_WIDTH           (CYCLES_WIDTH),
        .F_START_WIDTH          (3), 
        .TRUNCATED_ADDRESS_WIDTH(ADDRESS_WIDTH)
    ) u_phase_acc (
        .clk               (clk),
        .rst_n             (rst_n),
        .enable            (enable),      // <--- CONNECTED
        .FTW_start         (FTW_start),
        .cycles            (cycles),
        .FTW_step          (FTW_step),
        .truncated_address (phase_address),
        .valid_phase       (valid_phase)  // <--- CONNECTED
    );

    // --- 2. Quadrant Mapper ---
    first_quad_address #(
        .address_width(ADDRESS_WIDTH)
    ) u_quad_mapper (
        .truncated_address(phase_address),
        .first_address    (mapped_address),
        .negative_flag    (neg_flag),
        .clk              (clk)
    );

    // --- 3. Look-Up Table (Quarter Cycle) ---
    LUT #(
        .memory_width(MEMORY_WIDTH),
        .memory_depth(ADDRESS_WIDTH-2) 
    ) u_lut (
        .clk          (clk),
        .rst_n        (rst_n),
        .address      (mapped_address[13:0]), 
        .amplitude    (amplitude_out)
    );

    // --- 4. Negative MUX ---
    negative_mux mux (
        .amplitude_out   (amplitude_out),
        .negative_flag   (neg_flag),
        .final_amplitude (final_amplitude)
    );

    // --- 5. Valid Signal Pipeline Alignment ---
    // The LUT read and the negative_flag register both take 1 clock cycle.
    // We delay the valid signal by 1 clock cycle so it aligns perfectly with the final output.
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            valid_out_reg <= 1'b0;
        end else begin
            valid_out_reg <= valid_phase;
        end
    end
    
    assign valid_out = valid_out_reg;

endmodule

// `timescale 1ns / 1ps

// module dds_top #(
//     parameter TUNING_WORD_WIDTH = 32,
//     parameter CYCLES_WIDTH      = 13,
//     parameter ADDRESS_WIDTH     = 16, // Truncated phase width
//     parameter MEMORY_WIDTH      = 8   // DAC resolution
// )(
//     input  wire                         clk,
//     input  wire                         rst_n,
//     input  wire [TUNING_WORD_WIDTH-1:0]       FTW_start,
//     input  wire [CYCLES_WIDTH-1:0]      cycles,
//     input  wire [TUNING_WORD_WIDTH-1:0] FTW_step,
//     output wire [MEMORY_WIDTH-1:0]      final_amplitude
// );

//     // --- Internal Wires ---
//     wire [ADDRESS_WIDTH-1:0] phase_address;
//     wire [ADDRESS_WIDTH-1:0] mapped_address;
//     wire                     neg_flag;
//     wire [MEMORY_WIDTH-1:0]      amplitude_out;

//     // --- 1. Phase Accumulator ---
//     PHASE_ACC #(
//         .TUNING_WORD_WIDTH      (TUNING_WORD_WIDTH),
//         .CYCLES_WIDTH           (CYCLES_WIDTH),
//         .F_START_WIDTH          (3), 
//         .TRUNCATED_ADDRESS_WIDTH(ADDRESS_WIDTH)
//     ) u_phase_acc (
//         .FTW_start         (FTW_start),
//         .clk               (clk),
//         .rst_n             (rst_n),
//         .cycles            (cycles),
//         .FTW_step          (FTW_step),
//         .truncated_address(phase_address) 
//     );

//     // --- 2. Quadrant Mapper ---
//     first_quad_address #(
//         .address_width(ADDRESS_WIDTH)
//     ) u_quad_mapper (
//         .truncated_address(phase_address),
//         .first_address    (mapped_address),
//         .negative_flag    (neg_flag),
//         .clk(clk)
//     );

//     // --- 3. Look-Up Table (Quarter Cycle) ---
//     LUT #(
//         .memory_width(MEMORY_WIDTH),
//         .memory_depth(ADDRESS_WIDTH-2) 
//     ) u_lut (
//         .clk          (clk),
//         .rst_n        (rst_n),
//         .address      (mapped_address[13:0]), 
//         .amplitude    (amplitude_out)
//     );

//     // --- 4. Negative MUX ---
//     negative_mux mux (
//         .amplitude_out (amplitude_out),
//         .negative_flag (neg_flag),
//         .final_amplitude (final_amplitude)
//     );

// endmodule