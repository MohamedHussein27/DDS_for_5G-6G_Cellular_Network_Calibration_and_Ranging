// Analog Devices 
// GP Ain-shams University
// Butterfly unit for 4096-FFT 


/* Description ..........
    It uses MUX logic to correctly route the data into the delay feedback buffer or directly to the next stage based on the sel signal.
*/


module fft_stage #(
    parameter WL = 16,
    parameter DELAY_LEN = 2048,
    parameter ROM_DEPTH = 2048,
    parameter FILE_REAL = "",
    parameter FILE_IMAG = ""
)(
    input wire clk,
    input wire rst_n,
    input wire sel,
    input wire [(ROM_DEPTH == 1) ? 0 : $clog2(ROM_DEPTH)-1 : 0] addr,
    input wire signed [WL-1:0] in_real,
    input wire signed [WL-1:0] in_imag,
    output wire signed [WL-1:0] out_real,
    output wire signed [WL-1:0] out_imag
);

    wire signed [WL-1:0] delay_in_re, delay_in_im;
    wire signed [WL-1:0] delay_out_re, delay_out_im;
    wire signed [WL-1:0] bf_a_re, bf_a_im, bf_b_re, bf_b_im;
    wire signed [WL-1:0] mult_out_re, mult_out_im;
    wire signed [WL-1:0] twiddle_re, twiddle_im;

    // --- MUX Logic (Derived from standard SDF and our diagram) ---
    // MUX 1: Input to Delay Buffer (sel=0 takes x_in, sel=1 takes mult_out)
    assign delay_in_re = sel ? mult_out_re : in_real;
    assign delay_in_im = sel ? mult_out_im : in_imag;
    
    // MUX 2: Output of Stage (sel=0 takes delay_out, sel=1 takes Butterfly Add path)
    assign out_real = sel ? bf_a_re : delay_out_re;
    assign out_imag = sel ? bf_a_im : delay_out_im;

    // --- Sub-module Instantiations ---
    delayfeedback #(
        .WL(WL), 
        .L(DELAY_LEN) // Update your delayfeedback to use L directly
    ) delay_inst (
        .clk(clk),
        .rst_n(rst_n),
        .data_in_real(delay_in_re),
        .data_in_imag(delay_in_im),
        .data_out_real(delay_out_re),
        .data_out_imag(delay_out_im)
    );

    butterfly #(.WL(WL)) bf_inst (
        .in1_real(delay_out_re), // Butterfly 'a' gets delayed data
        .in1_imag(delay_out_im),
        .in2_real(in_real),      // Butterfly 'b' gets fresh input data
        .in2_imag(in_imag),
        .a_real(bf_a_re),        // Sum (goes to MUX 2)
        .a_imag(bf_a_im),
        .b_real(bf_b_re),        // Difference (goes to Multiplier)
        .b_imag(bf_b_im)
    );

    twiddle_rom #(
        .WL(WL),
        .DEPTH(ROM_DEPTH), 
        .FILE_REAL(FILE_REAL),
        .FILE_IMAG(FILE_IMAG)
    ) rom_inst (
        .addr_a('b0),                       // Reverted back to a simple, safe zero!
        .tw_re_a(),                         // Keep empty (Floating output)
        .tw_im_a(),                         // Keep empty (Floating output)
        .addr_b(addr),        
        .tw_re_b(twiddle_re), 
        .tw_im_b(twiddle_im)  
    );

    multiplier #(.WL(WL)) mult_inst (
        .re1(bf_b_re),
        .im1(bf_b_im),
        .re2(twiddle_re),
        .im2(twiddle_im),
        .re_out(mult_out_re),
        .im_out(mult_out_im)
    );

endmodule