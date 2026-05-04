`timescale 1ns / 1ps

module dds_config_regmap (
    input  wire        clk,
    input  wire        rst_n,

    // ==========================================
    // Simple Memory-Mapped Bus Interface
    // ==========================================
    input  wire [3:0]  addr,      // 4-bit address space
    input  wire        wr_en,     // Write Enable
    input  wire [31:0] wr_data,   // Data to write
    input  wire        rd_en,     // Read Enable
    output reg  [31:0] rd_data,   // Data read out

    // ==========================================
    // Direct Hardware Outputs (To DDS Module)
    // ==========================================
    output wire [31:0] FTW_start_out,
    output wire [31:0] FTW_step_out,
    output wire [12:0] cycles_out,
    
    // NEW: Signal to upper layer to enable the DDS
    output wire        dds_ready_flag 
);

    // ------------------------------------------
    // Address Memory Map Offsets
    // ------------------------------------------
    localparam ADDR_FTW_START = 4'h0; // Offset 0x0
    localparam ADDR_FTW_STEP  = 4'h4; // Offset 0x4
    localparam ADDR_CYCLES    = 4'h8; // Offset 0x8

    // ------------------------------------------
    // Internal Registers
    // ------------------------------------------
    reg [31:0] reg_ftw_start;
    reg [31:0] reg_ftw_step;
    reg [12:0] reg_cycles;

    // Continuously drive the configuration to the DDS
    assign FTW_start_out = reg_ftw_start;
    assign FTW_step_out  = reg_ftw_step;
    assign cycles_out    = reg_cycles;

    // ------------------------------------------
    // Write Logic & Reset Defaults
    // ------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            reg_ftw_start <= 32'd0;
            reg_ftw_step  <= 32'd426666;
            reg_cycles    <= 13'd4096;
        end else if (wr_en) begin
            case (addr)
                ADDR_FTW_START: reg_ftw_start <= wr_data;
                ADDR_FTW_STEP:  reg_ftw_step  <= wr_data;
                ADDR_CYCLES:    reg_cycles    <= wr_data[12:0];
            endcase
        end
    end

    // ------------------------------------------
    // Read Logic & Read Counter for Enable
    // ------------------------------------------
    reg [1:0] read_count;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            rd_data    <= 32'd0;
            read_count <= 2'd0;
        end else if (rd_en) begin
            // Track the number of reads up to 3
            if (read_count < 2'd3) begin
                read_count <= read_count + 1;
            end
            
            // Standard Read Logic
            case (addr)
                ADDR_FTW_START: rd_data <= reg_ftw_start;
                ADDR_FTW_STEP:  rd_data <= reg_ftw_step;
                ADDR_CYCLES:    rd_data <= {19'd0, reg_cycles}; 
                default:        rd_data <= 32'hDEADBEEF;       
            endcase
        end
    end

    // Assert ready flag once 3 reads have occurred
    assign dds_ready_flag = (read_count == 2'd3);

endmodule