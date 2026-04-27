`timescale 1ns / 1ps

module tb_ifft_2048();
    parameter WL = 16;
    parameter N = 2048;

    reg clk;
    reg rst_n;
    reg valid_in;
    reg signed [WL-1:0] in_real;
    reg signed [WL-1:0] in_imag;

    wire valid_out;
    wire signed [WL-1:0] out_real;
    wire signed [WL-1:0] out_imag;

    integer file_in_re, file_in_im, file_out_re, file_out_im;
    integer i, scan_re, scan_im;

    // Instantiate Top Module
    ifft_2048_top #(.WL(WL)) uut (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(valid_in),
        .in_real(in_real),
        .in_imag(in_imag),
        .valid_out(valid_out),
        .out_real(out_real),
        .out_imag(out_imag)
    );

    // Clock Generation
    initial begin
        clk = 0;
        forever #5 clk = ~clk;
    end

    // Dumps for Debugging
    initial begin
        $dumpfile("dump.vcd");
        $dumpvars(0, ifft_2048_top);
    end

    // Stimulus and Record
    initial begin
        // 1. Initialize
        rst_n = 0;
        valid_in = 0;
        in_real = 0;
        in_imag = 0;

        file_in_re = $fopen("test_data/input_real.txt", "r");
        file_in_im = $fopen("test_data/input_imag.txt", "r");
        file_out_re = $fopen("test_data/output_real.txt", "w");
        file_out_im = $fopen("test_data/output_imag.txt", "w");

        if (file_in_re == 0 || file_in_im == 0) begin
            $display("ERROR: Input files not found!");
            $finish;
        end

        // 2. Hardware Wake-up / Reset
        repeat(10) @(posedge clk);
        @(negedge clk);
        rst_n = 1;
        repeat(5) @(negedge clk);

        // 3. Feed the 2048 Data Samples
        for (i = 0; i < N; i = i + 1) begin
            valid_in = 1;
            scan_re = $fscanf(file_in_re, "%d\n", in_real);
            scan_im = $fscanf(file_in_im, "%d\n", in_imag);
            
            // Record output unconditionally (MATLAB will auto-align the delay)
            $fdisplay(file_out_re, "%d", out_real);
            $fdisplay(file_out_im, "%d", out_imag);
            
            @(negedge clk);
        end

        // 4. The Zero-Flush
        // Keep valid_in = 1, but feed zeros. This safely pushes the remaining 
        // math through the 2,047 delay registers without freezing the pipeline.
        in_real = 0;
        in_imag = 0;
        valid_in = 0; 
        
        for (i = 0; i < 2500; i = i + 1) begin
            $fdisplay(file_out_re, "%d", out_real);
            $fdisplay(file_out_im, "%d", out_imag);
            @(negedge clk);
        end

        // 5. Close and Finish
        $fclose(file_in_re);
        $fclose(file_in_im);
        $fclose(file_out_re);
        $fclose(file_out_im);
        $display("Simulation Complete. Output saved.");
        $finish;
    end
endmodule