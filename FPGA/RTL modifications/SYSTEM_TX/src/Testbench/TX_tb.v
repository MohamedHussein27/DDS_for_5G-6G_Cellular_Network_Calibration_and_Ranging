`timescale 1ns / 1ps

module tb_TX_wrapper();

    // 1. Signals declaration
    reg  clk;
    reg  rst_n;
    reg  dds_enable;
    
    // Only the flags and valid signal remain as outputs!
    wire tx_valid;
    wire real_match_flag;
    wire imag_match_flag;

    // 2. Instantiate the wrapper
    TX_wrapper #(
        .WL(16),
        .N(4096),
        .DDS_W(8)
    ) uut (
        .clk(clk),
        .rst_n(rst_n),
        .dds_enable(dds_enable),
        .tx_valid(tx_valid),
        .real_match_flag(real_match_flag),
        .imag_match_flag(imag_match_flag)
    );

    parameter CLK_PERIOD = 10;
    
    // 3. Clock Generation (100 MHz)
    initial begin
        clk = 0;
        forever #(CLK_PERIOD/2) clk = ~clk; 
    end

    // 4. Test Sequence
    initial begin
        // Initialize Inputs
        rst_n = 0;
        dds_enable = 0;

        // Wait 100 ns for global reset to finish
        #100;
        
        // Release reset
        rst_n = 1;
        #20;
        
        $display("[%0t] Starting DDS Configuration Phase...", $time);
        
        // Assert dds_enable to start the counter and config sequence
        dds_enable = 1;
        
        // The dds_counter takes a few cycles to step through 0, 1, 2, 3
        // We hold enable high to let it count and lock the configuration
        #(CLK_PERIOD * 10); 
        dds_enable = 0;
        
        $display("[%0t] Config loaded. Awaiting valid TX output...", $time);
        $display("[%0t] Note: This will take >12,000 clock cycles due to FFT + Bit Reversal + IFFT pipeline depth.", $time);
        
        // Wait for the pipeline to finish and assert valid
        wait(tx_valid == 1'b1);
        $display("[%0t] SUCCESS: tx_valid asserted! Pipeline is outputting data.", $time);
        
        // Let it run for the full 4096-point frame to allow comparators to calculate majority vote
        #(CLK_PERIOD * 25000); 
        
        // Check final results
        $display("========================================");
        if (real_match_flag)
            $display("[%0t] REAL DATA MATCH: PASSED", $time);
        else
            $display("[%0t] REAL DATA MATCH: FAILED", $time);
            
        if (imag_match_flag)
            $display("[%0t] IMAG DATA MATCH: PASSED", $time);
        else
            $display("[%0t] IMAG DATA MATCH: FAILED", $time);
        $display("========================================");
        
        $display("[%0t] Simulation finished.", $time);
        $stop;
    end

endmodule