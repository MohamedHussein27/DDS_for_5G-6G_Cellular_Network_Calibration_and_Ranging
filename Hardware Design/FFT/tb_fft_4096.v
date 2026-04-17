`timescale 1ns / 1ps

module tb_fft_4096();

    // Parameters
    parameter WL = 16;
    parameter N = 4096;
    
    // Clock and Reset
    reg clk;
    reg rst_n;
    
    // Inputs
    reg signed [WL-1:0] in_real;
    reg signed [WL-1:0] in_imag;
    
    // Outputs
    wire signed [WL-1:0] out_real;
    wire signed [WL-1:0] out_imag;
    
    // File Pointers
    integer file_in_re, file_in_im;
    integer file_out_re, file_out_im;
    integer scan_re, scan_im;
    
    // Cycle Tracker
    integer cycle_count;

    // Instantiate the Unit Under Test (UUT)
    fft_4096_top #(.WL(WL)) uut (
        .clk(clk),
        .rst_n(rst_n),
        .in_real(in_real),
        .in_imag(in_imag),
        .out_real(out_real),
        .out_imag(out_imag)
    );

    // Clock Generation (100 MHz)
    initial begin
        clk = 0;
        forever #5 clk = ~clk; 
    end

    // Test Sequence
    initial begin
        // 1. Initialize Inputs
        in_real = 0;
        in_imag = 0;
        cycle_count = 0;
        
        // 2. Open Files
        // Make sure these files exist in your simulation directory!
        file_in_re = $fopen("test_data/input_real.txt", "r");
        file_in_im = $fopen("test_data/input_imag.txt", "r");
        
        file_out_re = $fopen("test_data/output_real.txt", "w");
        file_out_im = $fopen("test_data/output_imag.txt", "w");

        if (file_in_re == 0 || file_in_im == 0) begin
            $display("ERROR: Could not open input files.");
            $finish;
        end

        // 3. Apply Reset
        rst_n = 0;
        #20; 
        rst_n = 1;
        
        // 4. Feed Data into the FFT (Two-Phase Method)
        
        // Sync to the first falling edge to avoid all race conditions
        @(negedge clk); 
        
        // =========================================================
        // PHASE 1: Feed Exactly 4096 Samples
        // =========================================================
        for (cycle_count = 0; cycle_count < N; cycle_count = cycle_count + 1) begin
            // 1. Put data on the wire
            scan_re = $fscanf(file_in_re, "%d\n", in_real);
            scan_im = $fscanf(file_in_im, "%d\n", in_imag);
            
            // 2. Capture valid outputs (First valid output appears at N-1)
            if (cycle_count >= (N - 1)) begin
                $fdisplay(file_out_re, "%d", out_real);
                $fdisplay(file_out_im, "%d", out_imag);
            end
            
            // 3. Wait for next clock cycle
            @(negedge clk); 
        end
        
        // =========================================================
        // PHASE 2: Flush Pipeline (Feed Zeros)
        // =========================================================
        // Instantly force wires to 0 exactly after the 4096th sample
        in_real = 0;
        in_imag = 0;
        
        for (cycle_count = N; cycle_count <= (N * 2) + 20; cycle_count = cycle_count + 1) begin
            // 1. Continue capturing the remaining outputs
            if (cycle_count <= ((N * 2) - 2)) begin
                $fdisplay(file_out_re, "%d", out_real);
                $fdisplay(file_out_im, "%d", out_imag);
            end
            
            // 2. Wait for next clock cycle
            @(negedge clk);
        end
        
        // Let the simulation breathe before finishing
        #50;

        // 5. Clean up and end simulation
        $fclose(file_in_re);
        $fclose(file_in_im);
        $fclose(file_out_re);
        $fclose(file_out_im);
        
        $display("Simulation Complete. Check output files.");
        $stop;
    end

endmodule