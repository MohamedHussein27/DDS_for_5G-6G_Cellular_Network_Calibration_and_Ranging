module multiplier_algorithm #(parameter WL=16) (
    input wire clk,
    input wire rst_n,
    input wire signed [WL-1:0] re1, im1, re2, im2,
    output reg signed [WL-1:0] re_out, im_out
);
    // 32-bit wires for full precision
    wire signed [31:0] mul1, mul2, mul3, mul4;
    // 33-bit wires for Add/Sub overflow protection
    wire signed [32:0] sub_re, add_im; 
    // Q8.8 * Q2.14 = Q10.22

    multiplier_fft u_mult_rxr (
        .A (re1),
        .B (re2),
        .P (mul1)
    );

    multiplier_fft u_mult_ximxim (
        .A (im1), 
        .B (im2),
        .P (mul2)
    ); 

    multiplier_fft u_imxr (
        .A (im1),
        .B (re2),
        .P (mul3)
    );

    multiplier_fft u_rxim (
        .A (re1), 
        .B (im2), 
        .P (mul4)
    );

    assign sub_re = mul1 - mul2;
    assign add_im = mul3 + mul4;


    // Shift by 14 to return to Q8.8. 
    // Adding bit [13] provides nearest-integer rounding.
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            re_out <= 0;
            im_out <= 0;
        end else begin
            re_out <= sub_re[29:14] + sub_re[13]; 
            im_out <= add_im[29:14] + add_im[13];
        end
    end
endmodule