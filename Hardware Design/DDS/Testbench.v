`timescale 1ns / 1ps

module tb_dds_top;

    // --- System Signals ---
    reg         clk;
    reg         rst_n;
    reg   [31:0] FTW_start;
    reg   [12:0] cycles;
    reg   [31:0] FTW_step;
    wire  [7:0] final_amplitude; 

    integer file_id;

    // --- Clock Generation ---
    initial begin
        clk = 0;
        forever #5 clk = ~clk; // 10ns period -> 100 MHz clock frequency
    end

    // --- DUT Instantiation ---
    dds_top u_dut (
        .clk          (clk),
        .rst_n        (rst_n),
        .FTW_start    (FTW_start),
        .cycles       (cycles),
        .FTW_step     (FTW_step),
        .final_amplitude(final_amplitude)
    );

    // --- Stimulus and File I/O ---
    initial begin
        // Open file to save RTL output for MATLAB spectral analysis
        file_id = $fopen("rtl_output.txt", "w");
        
        // Initialize Inputs
        rst_n     = 0;
        cycles    = 13'd4096;      // Will run for 4096 clock cycles before resetting
        FTW_step  = 32'd426666;         // Set to 0 for a constant single-tone frequency (no chirping)
        
        // FTW determines output frequency: Fout = (FTW * Fclk) / 2^32
        // For 100MHz clock, FTW=174762666 yields approx 4.068 MHz
        FTW_start = 32'd0; 
        
        #20;
        rst_n = 1; // Release reset
        
        // Run simulation for 4100 cycles to capture a full data set
        repeat (4100) begin
            @(posedge clk);
            // Wait a fraction of a cycle to ensure data is stable before sampling
            #1; 
            
            // Cast to signed to ensure MATLAB reads negative values correctly, 
            // rather than large positive integers
            $fdisplay(file_id, "%d", $signed(final_amplitude)); 
        end

        $fclose(file_id);
        $display("Simulation complete. Output written to rtl_output.txt");
        $finish;
    end

endmodule