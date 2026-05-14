`timescale 1ns / 1ps

module RX_wrapper_tb();

    // =========================================================
    // PARAMETERS & SIGNALS
    // =========================================================
    parameter WL = 16;
    parameter N  = 4096;

    reg clk;
    reg rst_n;

    // Channel Input Stimulus
    reg                 rx_valid_in;

    // Channel outputs (for monitoring)
    wire real_match_flag_ofdm;
    wire imag_match_flag_ofdm;
    wire real_match_flag_radar;
    wire imag_match_flag_radar;



    // =========================================================
    // DUT INSTANTIATION (Corrected to RX_wrapper)
    // =========================================================
    RX_wrapper #(
        .WL(WL), .N(N)
    ) uut (
        .clk(clk),
        .rst_n(rst_n),
        .rx_valid_in(rx_valid_in),
        .real_match_flag_ofdm(real_match_flag_ofdm),
        .imag_match_flag_ofdm(imag_match_flag_ofdm),
        .real_match_flag_radar(real_match_flag_radar),
        .imag_match_flag_radar(imag_match_flag_radar)
    );

    // =========================================================
    // CLOCK GENERATION
    // =========================================================
    parameter CLK_PERIOD = 10; // 100 MHz
    
    initial begin
        clk = 0;
        forever #(CLK_PERIOD/2) clk = ~clk; // 10ns period (100 MHz)
    end

    // =========================================================
    // MAIN STIMULUS (50,000 Cycle Assert)
    // =========================================================
    initial begin
        // 1. Initialize Inputs
        rst_n       = 0;
        rx_valid_in = 0;

        // 2. Apply Reset
        #20;
        rst_n = 1;
        #25;

        // 3. Assert the enables
        $display("[%0t] Asserting enable signals...", $time);


        rx_valid_in = 1'b1;

        // 4. Wait for exactly 50,000 clock cycles
        //    (Dummy data added so the FFT is processing changing numbers)
        #(4096 * CLK_PERIOD); // Wait for 4097 cycles to ensure we cover the entire input sequence (0 to 4095) and one extra cycle for processing

        rx_valid_in = 1'b0;

        #(20000 * CLK_PERIOD); // Wait additional cycles to allow pipeline to flush and outputs to stabilize

        // 5. Deassert the enables and Stop
        $display("[%0t] 50,000 cycles reached. Deasserting enables.", $time);


        // 6. Give the pipeline a few extra cycles to flush out final data
        #200; 
        
        $display("[%0t] Simulation Complete.", $time);
        $stop;
    end

endmodule