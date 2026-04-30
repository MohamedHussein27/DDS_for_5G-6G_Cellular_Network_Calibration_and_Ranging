module wrapper #(
    parameter TUNING_WORD_WIDTH = 32,
    parameter CYCLES_WIDTH      = 13,
    parameter TRUNCATED_ADDRESS_WIDTH     = 16, // Truncated phase width
    parameter MEMORY_WIDTH      = 8,   // DAC resolution
    parameter COUNTER_WIDTH = 12,      // Fixed semicolon and added comma
    parameter NUM_SYMBOLS = 4096
) (
    input  wire                         clk,
    input  wire                         rst_n,
    input  wire                         enable,     
    output wire                         comparison_result           
);          

// Changed to localparam and added widths to match port connections
    localparam [TUNING_WORD_WIDTH-1:0] FTW_start = 32'd0;
    localparam [TUNING_WORD_WIDTH-1:0] FTW_step  = 32'd426666;
    localparam [CYCLES_WIDTH-1:0]      cycles    = 13'd4096;


    wire valid_out;
    wire [MEMORY_WIDTH-1:0] final_amplitude;
    wire [COUNTER_WIDTH-1:0] write_addr;
    wire signed [MEMORY_WIDTH-1:0] reference_out;

    dds_top #(
        .TUNING_WORD_WIDTH(TUNING_WORD_WIDTH),
        .CYCLES_WIDTH(CYCLES_WIDTH), 
        .ADDRESS_WIDTH(TRUNCATED_ADDRESS_WIDTH), // Name matched to dds_top
        .MEMORY_WIDTH(MEMORY_WIDTH)
    ) u_dds (
        .clk (clk),
        .rst_n (rst_n),
        .enable (enable),     
        .FTW_start (FTW_start),
        .cycles (cycles),
        .FTW_step (FTW_step),
        .valid_out (valid_out), // Added missing comma
        .final_amplitude (final_amplitude)
    );

    bram_write_counter #(
        .ADDR_WIDTH(COUNTER_WIDTH) // Parameter name matched to counter
    ) u_count (
        .clk (clk),
        .rst_n (rst_n),
        .valid_in (valid_out),
        .write_addr (write_addr)
    );

    Reference_output_mem u_ref_mem (
      .clka(clk),
      .addra(write_addr),
      .douta(reference_out) 
    );

    comparator #(.INPUT_WIDTH(MEMORY_WIDTH), .NUM_SYMBOLS(NUM_SYMBOLS), .SYMBOL(COUNTER_WIDTH)) u_comp (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(enable),
        .dds_in (final_amplitude),
        .ref_in (reference_out),
        .symbol_index (write_addr),
        .comparison_result (comparison_result)
    );


endmodule