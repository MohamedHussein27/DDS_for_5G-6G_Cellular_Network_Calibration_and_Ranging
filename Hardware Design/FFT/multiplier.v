// Analog Devices 
// GP Ain-shams University
// Multiplier block for 4096-FFT 


/* Description ..........
    multiplier take 2 inputs which are the twiddle factor from the twiddle ROM and the "lower" subtraction output of the butterfly unit
    and each one of them has a real and imaginary part, so we multiply 2 complex numbers so we need 4 multipliers and 2 adders to do this
    (re1 + jim1) * (re2 + jim2) = (re1*re2 - im1*im2) + j(im1*re2 + re1*im2)
*/


module multiplier #(parameter WL=16) (
    input wire signed [WL-1:0] re1, im1, re2, im2,
    output wire signed [WL-1:0] re_out, im_out
);
    // 32-bit wires for full precision
    wire signed [31:0] mul1, mul2, mul3, mul4;
    // 33-bit wires for Add/Sub overflow protection
    wire signed [32:0] sub_re, add_im; 

    // Q2.14 * Q2.14 = Q4.28
    assign mul1 = re1 * re2;
    assign mul2 = im1 * im2;
    assign mul3 = im1 * re2;
    assign mul4 = re1 * im2;

    assign sub_re = mul1 - mul2;
    assign add_im = mul3 + mul4;

    // Shift by 14 to return to Q2.14. 
    // Adding bit [13] provides nearest-integer rounding.
    assign re_out = sub_re[29:14] + sub_re[13]; 
    assign im_out = add_im[29:14] + add_im[13];

endmodule