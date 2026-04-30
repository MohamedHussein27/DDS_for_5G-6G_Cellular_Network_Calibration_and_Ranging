// Analog Devices 
// GP Ain-shams University
// Butterfly unit for 4096-FFT 


/* Description ..........
    It uses MUX logic to correctly route the data into the delay feedback buffer or directly to the next stage based on the sel signal.
*/


// Analog Devices 
// GP Ain-shams University
// FFT Stage for 2048-point FFT 

/* Description ..........
    It uses MUX logic to correctly route the data into the delay feedback buffer or directly to the next stage based on the sel signal.
    The valid_in signal acts as an enable for the delay feedback buffer to keep the pipeline synchronized.
*/

module ifft_stage_2048 #(
    parameter WL = 16,
    parameter DELAY_LEN = 1024, // Note: For a 2048 FFT, max delay is 1024
    parameter ROM_DEPTH = 1024
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in, // <--- UPDATED: Pipeline Enable / Valid signal
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

    // --- MUX Logic (Derived from standard SDF architecture) ---
    // MUX 1: Input to Delay Buffer (sel=0 takes fresh input, sel=1 takes multiplier result)
    assign delay_in_re = sel ? mult_out_re : in_real;
    assign delay_in_im = sel ? mult_out_im : in_imag;
    
    // MUX 2: Output of Stage (sel=0 takes delay_out, sel=1 takes Butterfly Add result)
    assign out_real = sel ? bf_a_re : delay_out_re;
    assign out_imag = sel ? bf_a_im : delay_out_im;

    // --- Sub-module Instantiations ---
    
    // Delay Feedback Buffer
    delayfeedback #(
        .WL(WL), 
        .L(DELAY_LEN) 
    ) delay_inst (
        .clk(clk),
        .rst_n(rst_n),
        .en(valid_in), // <--- UPDATED: Pass valid_in as enable
        .data_in_real(delay_in_re),
        .data_in_imag(delay_in_im),
        .data_out_real(delay_out_re),
        .data_out_imag(delay_out_im)
    );

    // Butterfly Unit
    butterfly #(.WL(WL)) bf_inst (
        .in1_real(delay_out_re), // Butterfly 'a' gets delayed data
        .in1_imag(delay_out_im),
        .in2_real(in_real),      // Butterfly 'b' gets fresh input data
        .in2_imag(in_imag),
        .a_real(bf_a_re),        // Sum
        .a_imag(bf_a_im),
        .b_real(bf_b_re),        // Difference
        .b_imag(bf_b_im)
    );

    // Twiddle Factor ROM for 2048-point FFT
    twiddlerom_2048 #(.WL(WL), .DEPTH(ROM_DEPTH)) rom_inst (
        .addr_a(1'b0), // for IFFT we use the conjugate of the twiddle so we used another address to read the conjugate value from the ROM
        .addr_b(addr),
        .W_real_a(),  .W_img_a(),
        .W_real_b(twiddle_re), .W_img_b(twiddle_im)
    );

    // Complex Multiplier
    multiplier #(.WL(WL)) mult_inst (
        .re1(bf_b_re),
        .im1(bf_b_im),
        .re2(twiddle_re),
        .im2(twiddle_im),
        .re_out(mult_out_re),
        .im_out(mult_out_im)
    );

endmodule