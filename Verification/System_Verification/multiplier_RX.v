`timescale 1ns / 1ps

module multiplier_RX #(parameter WL=16) (
    input wire signed [WL-1:0] re1, im1, re2, im2,
    output wire signed [WL-1:0] re_out, im_out
);
    // 32-bit wires for full precision
    wire signed [31:0] mul1, mul2, mul3, mul4;
    
    // 33-bit wires for Add/Sub overflow protection
    wire signed [32:0] sub_re, add_im; 

    // Q8.8 * Q2.14 = Q10.22
    assign mul1 = re1 * re2;
    assign mul2 = im1 * im2;
    assign mul3 = im1 * re2;
    assign mul4 = re1 * im2;

    // Addition requires 1 extra bit -> Q11.22
    assign sub_re = mul1 - mul2;
    assign add_im = mul3 + mul4;

    // Shift right by 17 to return to Q11.5. 
    // Bits [32:17] map exactly to the 16-bit output.
    // Adding bit [16] provides nearest-integer rounding.
    assign re_out = sub_re[32:17] + sub_re[16]; 
    assign im_out = add_im[32:17] + add_im[16];

endmodule