// Analog Devices 
// GP Ain-shams University
// Butterfly unit for 4096-FFT 


/* Description ..........
    It uses MUX logic to correctly route the data into the delay feedback buffer or directly to the next stage based on the sel signal.
*/



module fft_stage #(
    parameter WL = 16,
    parameter DELAY_LEN = 2048,
    parameter ROM_DEPTH = 2048
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in, // <--- NEW: Pipeline Enable
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

    assign delay_in_re = sel ? mult_out_re : in_real;
    assign delay_in_im = sel ? mult_out_im : in_imag;

    assign out_real = sel ? bf_a_re : delay_out_re;
    assign out_imag = sel ? bf_a_im : delay_out_im;

    delayfeedback #(
        .WL(WL), 
        .L(DELAY_LEN) 
    ) delay_inst (
        .clk(clk),
        .rst_n(rst_n),
        .en(valid_in), // <--- NEW: Pass enable to buffer
        .data_in_real(delay_in_re),
        .data_in_imag(delay_in_im),
        .data_out_real(delay_out_re),
        .data_out_imag(delay_out_im)
    );

    butterfly #(.WL(WL)) bf_inst (
        .in1_real(delay_out_re), .in1_imag(delay_out_im),
        .in2_real(in_real),      .in2_imag(in_imag),
        .a_real(bf_a_re),        .a_imag(bf_a_im),
        .b_real(bf_b_re),        .b_imag(bf_b_im)
    );

    twiddlerom #(.WL(WL), .DEPTH(ROM_DEPTH)) rom_inst (
        .addr_a(0),
        .addr_b(addr),
        .W_real_a(), 
        .W_img_a(),
        .W_real_b(twiddle_re),           
        .W_img_b(twiddle_im)
    );

    multiplier #(.WL(WL)) mult_inst (
        .re1(bf_b_re),    .im1(bf_b_im),
        .re2(twiddle_re), .im2(twiddle_im),
        .re_out(mult_out_re), .im_out(mult_out_im)
    );
endmodule
