`timescale 1ns / 1ps

module tb_FFT_wrapper;

    // Parameters
    parameter WL = 16;
    parameter ADDR_WIDTH = 12;
    parameter INPUT_WIDTH = 16;
    parameter NUM_SYMBOLS = 4096;
    parameter SYMBOL = 12;
    
    parameter CLK_PERIOD = 2; // Assuming 500 MHz based on #1 toggles

    // Inputs
    reg clk;
    reg rst_n;
    reg start;

    // Outputs
    wire real_comparison_result;
    wire imag_comparison_result;

    // Instantiate the Unit Under Test (UUT)
    FFT_wrapper #(
        .WL(WL),
        .ADDR_WIDTH(ADDR_WIDTH),
        .INPUT_WIDTH(INPUT_WIDTH),
        .NUM_SYMBOLS(NUM_SYMBOLS),
        .SYMBOL(SYMBOL)
    ) uut (
        .clk(clk),
        .rst_n(rst_n),
        .start(start),
        .real_comparison_result(real_comparison_result),
        .imag_comparison_result(imag_comparison_result)
    );

    // Clock generation
    initial begin
        clk = 0;
        forever #(CLK_PERIOD / 2.0) clk = ~clk;
    end

    // Test Stimulus
    initial begin
        // 1. Initialize Inputs
        rst_n = 0;
        start = 0;

        // 2. Wait for global reset to finish
        #(CLK_PERIOD * 50);
        
        // 3. Release active-low reset
        rst_n = 1;
        
        // Wait a few clock cycles before starting
        #(CLK_PERIOD * 5);
        
        // 4. Assert start signal to begin streaming 4096 samples
        $display("[%0t] Starting FFT Data Stream (4096 cycles)...", $time);
        start = 1;
        
        // Keep start high for exactly 4096 clock cycles
        #(CLK_PERIOD * 4096);
        start = 0;
        
        $display("[%0t] Input stream finished. Waiting for pipeline to auto-flush...", $time);
        
        // 5. Wait for computation to complete.
        // The new control unit takes 4107 cycles to flush. We wait 4500 to be safe.
        #(CLK_PERIOD * 4500);
        
        // 6. End simulation and print final results
        if (error_count == 0) begin
            $display("---------------------------------------------------------");
            $display(">>> TEST PASSED: All 4096 outputs matched the reference! <<<");
            $display("---------------------------------------------------------");
        end else begin
            $display("---------------------------------------------------------");
            $display(">>> TEST FAILED: Found %0d mismatches! <<<", error_count);
            $display("---------------------------------------------------------");
        end
        $finish;
    end
    
    // ---------------------------------------------------------
    // Error Monitoring
    // ---------------------------------------------------------
    integer error_count = 0;
    
    always @(posedge clk) begin
        // Only check the comparison flags when uut.valid_out is HIGH!
        if (rst_n && uut.valid_out) begin
            // Assuming 0 means mismatch in your comparator
            if (real_comparison_result === 1'b0 || imag_comparison_result === 1'b0) begin
                $display("Mismatch detected at ROM Index %0d! Real_Match=%b, Imag_Match=%b", 
                         uut.write_addr, real_comparison_result, imag_comparison_result);
                error_count = error_count + 1;
            end
        end
    end

endmodule