`timescale 1ns / 1ps

module dds_top #(
    parameter TUNING_WORD_WIDTH = 32, // Width of the phase accumulator
    parameter CYCLES_WIDTH      = 13, // Counter width for duration control
    parameter ADDRESS_WIDTH     = 16, // Phase bits kept after truncation (determines phase resolution)
    parameter MEMORY_WIDTH      = 8   // Width of the LUT data (determines DAC amplitude resolution)
)(
    input  wire                         clk,
    input  wire                         rst_n,
    input  wire [TUNING_WORD_WIDTH-1:0] FTW_start, // Initial Frequency Tuning Word
    input  wire [CYCLES_WIDTH-1:0]      cycles,    // Number of cycles before resetting/wrapping
    input  wire [TUNING_WORD_WIDTH-1:0] FTW_step,  // Value added to FTW every clock cycle (for chirping/FM)
    output wire [MEMORY_WIDTH-1:0]      amplitude_out
);

    // --- Internal Wires ---
    // Carries the truncated phase from the accumulator
    wire [ADDRESS_WIDTH-1:0] phase_address;
    
    // Carries the phase folded into the 1st quadrant (0 to pi/2)
    wire [ADDRESS_WIDTH-1:0] mapped_address;
    
    // Flag indicating if we are in the 3rd or 4th quadrant (requires amplitude inversion)
    wire                     neg_flag;

    // --- 1. Phase Accumulator ---
    // Generates a digital ramp representing the instantaneous phase
    PHASE_ACC #(
        .TUNING_WORD_WIDTH      (TUNING_WORD_WIDTH),
        .CYCLES_WIDTH           (CYCLES_WIDTH),
        .F_START_WIDTH          (3), 
        .TRUNCATED_ADDRESS_WIDTH(ADDRESS_WIDTH)
    ) u_phase_acc (
        .clk               (clk),
        .rst_n             (rst_n),
        .FTW_start         (FTW_start),
        .cycles            (cycles),
        .FTW_step          (FTW_step),
        .truncated_address (phase_address) 
    );

    // --- 2. Quadrant Mapper ---
    // Exploits sine wave symmetry to reduce LUT memory size by a factor of 4
    first_quad_address #(
        .address_width(ADDRESS_WIDTH)
    ) u_quad_mapper (
        .truncated_address(phase_address),
        .first_address    (mapped_address),
        .negative_flag    (neg_flag)
    );

    // --- 3. Look-Up Table (Quarter Cycle) ---
    // Stores only 0 to 90 degrees of the sine wave
    LUT #(
        .memory_width(MEMORY_WIDTH),
        .memory_depth(ADDRESS_WIDTH-2) // Depth is 1/4 of the total phase space
    ) u_lut (
        .clk          (clk),
        .rst_n        (rst_n),
        .address      (mapped_address[13:0]), // Use lower bits for quarter-wave address mapping
        .negative_flag(neg_flag),
        .amplitude    (amplitude_out)
    );

endmodule