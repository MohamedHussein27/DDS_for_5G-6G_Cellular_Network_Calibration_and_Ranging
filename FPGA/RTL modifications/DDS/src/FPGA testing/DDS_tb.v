`timescale 1ns / 1ps

module tb_wrapper();

    // --- Parameters ---
    parameter TUNING_WORD_WIDTH = 32;
    parameter CYCLES_WIDTH      = 13;
    parameter TRUNCATED_ADDRESS_WIDTH = 16;
    parameter MEMORY_WIDTH      = 8;
    parameter COUNTER_WIDTH     = 12;
    parameter NUM_SYMBOLS       = 4096;
    
    // --- Clock Period ---
    localparam CLK_PERIOD = 10;

    // --- Signals ---
    reg  clk;
    reg  rst_n;
    reg  enable;
    
    wire comparison_result; 
    
    // --- Instantiation ---
    wrapper #(
        .TUNING_WORD_WIDTH(TUNING_WORD_WIDTH),
        .CYCLES_WIDTH(CYCLES_WIDTH),
        .TRUNCATED_ADDRESS_WIDTH(TRUNCATED_ADDRESS_WIDTH),
        .MEMORY_WIDTH(MEMORY_WIDTH),
        .COUNTER_WIDTH(COUNTER_WIDTH),
        .NUM_SYMBOLS(NUM_SYMBOLS)
    ) uut (
        .clk(clk),
        .rst_n(rst_n),
        .enable(enable),
        .comparison_result(comparison_result)
    );

    // --- Clock Generation ---
    initial begin
        clk = 0;
        forever #(CLK_PERIOD/2) clk = ~clk;
    end

    // --- Stimulus Generation ---
    initial begin
        rst_n  = 1'b0; 
        enable = 1'b0;
        
        #100;
        
        rst_n = 1'b1;
        
        #(CLK_PERIOD * 5);
        
        $display("Starting DDS Chirp...");
        enable = 1'b1;
        
        // FIXED: Using NUM_SYMBOLS instead of the deleted 'cycles' variable
        #(CLK_PERIOD * (2 * NUM_SYMBOLS));
        
        enable = 1'b0;
        $display("Chirp Complete.");
        
        #50;
        if (comparison_result == 1'b1) begin
            $display("========================================");
            $display("             TEST PASSED!               ");
            $display("   No mismatches found in the stream.   ");
            $display("========================================");
        end else begin
            $display("========================================");
            $display("             TEST FAILED!               ");
            $display("  Mismatch detected between DDS and ROM.");
            $display("========================================");
        end
        
        $finish;
    end

endmodule