`timescale 1ns / 1ps

module tb_fft_4096();

    // Parameters
    parameter WL = 16;
    parameter N = 4096;
    
    // Clock and Reset
    reg clk;
    reg rst_n;
    
    // Handshaking
    reg  valid_in;
    wire valid_out;
    
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
    integer i;

    // Instantiate the Unit Under Test (UUT)
    fft_4096_top #(.WL(WL)) uut (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(valid_in),   // New port
        .in_real(in_real),
        .in_imag(in_imag),
        .valid_out(valid_out), // New port
        .out_real(out_real),
        .out_imag(out_imag)
    );

    // Clock Generation (100 MHz)
    initial begin
        clk = 0;
        forever #5 clk = ~clk; 
    end

    // Monitor: Automatically capture output whenever valid_out is high
    initial begin
        wait(rst_n);
        forever begin
            @(negedge clk);
            if (valid_out) begin
                $fdisplay(file_out_re, "%d", out_real);
                $fdisplay(file_out_im, "%d", out_imag);
            end
        end
    end

    // Test Sequence
    initial begin
        // 1. Initialize Inputs
        in_real = 0;
        in_imag = 0;
        valid_in = 0;
        
        // 2. Open Files
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
        #50; 
        rst_n = 1;
        
        // Sync to falling edge
        @(negedge clk); 
        
        // =========================================================
        // PHASE 1: Feed Exactly 4096 Samples
        // =========================================================
        $display("Starting Data Ingestion...");
        valid_in = 1;
        for (i = 0; i < N; i = i + 1) begin
            scan_re = $fscanf(file_in_re, "%d\n", in_real);
            scan_im = $fscanf(file_in_im, "%d\n", in_imag);
            @(negedge clk); 
        end
        
        // =========================================================
        // PHASE 2: End Ingestion and Wait for Auto-Flush
        // =========================================================
        valid_in = 0;
        in_real  = 0;
        in_imag  = 0;
        $display("Ingestion Complete. Waiting for Auto-Flush to finish...");

        // The simulation will wait until valid_out goes low again
        // this signifies that the internal state machine has finished the frame.
        wait(valid_out == 1); // Wait for output to start
        wait(valid_out == 0); // Wait for output to finish

        // Let the simulation breathe
        #100;

        // 5. Clean up
        $fclose(file_in_re);
        $fclose(file_in_im);
        $fclose(file_out_re);
        $fclose(file_out_im);
        
        $display("Simulation Complete. Output files generated.");
        $stop;
    end

endmodule
