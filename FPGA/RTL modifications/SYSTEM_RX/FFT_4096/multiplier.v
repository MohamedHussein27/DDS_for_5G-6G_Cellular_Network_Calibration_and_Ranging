module multiplier #(parameter WL=16) (
    input wire clk,
    input wire rst_n,
    input wire en,
    input wire signed [WL-1:0] re1, im1, re2, im2,
    output reg signed [WL-1:0] re_out, im_out
);
    wire signed [31:0] mul1, mul2, mul3, mul4;
    wire signed [32:0] sub_re, add_im; 

    // assign mul1 = re1 * re2;
    // assign mul2 = im1 * im2;
    // assign mul3 = im1 * re2;
    // assign mul4 = re1 * im2;

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

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            re_out <= 0;
            im_out <= 0;
        end else if (en) begin
            re_out <= sub_re[29:14] + sub_re[13]; 
            im_out <= add_im[29:14] + add_im[13];
        end
    end
endmodule

