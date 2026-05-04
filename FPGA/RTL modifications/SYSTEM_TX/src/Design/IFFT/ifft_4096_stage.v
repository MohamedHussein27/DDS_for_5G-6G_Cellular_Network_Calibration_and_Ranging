module ifft_stage_4096 #(
    parameter WL = 16,
    parameter DELAY_LEN = 2048,
    parameter ROM_DEPTH = 2048
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in, 
    input wire sel,
//    input wire [(ROM_DEPTH == 1) ? 0 : $clog2(ROM_DEPTH)-1 : 0] addr,
    input wire signed [WL-1:0] in_real,
    input wire signed [WL-1:0] in_imag,
    input wire signed [WL-1:0] twiddle_re,
    input wire signed [WL-1:0] twiddle_im,
    output wire signed [WL-1:0] out_real,
    output wire signed [WL-1:0] out_imag
);

//    localparam MASTER_DEPTH = 2048;
//    localparam STRIDE = MASTER_DEPTH / ROM_DEPTH;

    wire signed [WL-1:0] delay_in_re, delay_in_im;
    wire signed [WL-1:0] delay_out_re, delay_out_im;
    wire signed [WL-1:0] bf_a_re, bf_a_im, bf_b_re, bf_b_im;
    wire signed [WL-1:0] mult_out_re, mult_out_im;
    wire signed [WL-1:0] twiddle_im_mult;

    assign twiddle_im_mult = -twiddle_im;

    delayfeedback #(.WL(WL), .L(DELAY_LEN - 1)) delay_inst (
        .clk(clk), .rst_n(rst_n), .en(valid_in), 
        .data_in_real(delay_in_re), .data_in_imag(delay_in_im),
        .data_out_real(delay_out_re), .data_out_imag(delay_out_im)
    );

    butterfly #(.WL(WL)) bf_inst (
        .in1_real(delay_out_re), .in1_imag(delay_out_im),
        .in2_real(in_real),      .in2_imag(in_imag),
        .a_real(bf_a_re),        .a_imag(bf_a_im),
        .b_real(bf_b_re),        .b_imag(bf_b_im)
    );

    // twiddle_real u_twiddle_real (
    //     .addra(addr * STRIDE), .clka(clk),
    //     .douta(twiddle_re)
    // );

    // twiddle_imag u_twiddle_imag (
    //     .addra(addr * STRIDE), .clka(clk),
    //     .douta(twiddle_im)
    // );

    // twiddlerom_4096 #(.WL(WL), .DEPTH(ROM_DEPTH)) rom_inst (
    //     .addr(addr), .clk(clk),
    //     .W_real(twiddle_re), .W_img(twiddle_im)
    // );

    multiplier #(.WL(WL)) mult_inst (
        .clk(clk), .rst_n(rst_n), .en(valid_in),
        .re1(bf_b_re),    .im1(bf_b_im),
        .re2(twiddle_re), .im2(twiddle_im_mult),
        .re_out(mult_out_re), .im_out(mult_out_im)
    );

    reg signed [WL-1:0] bf_a_re_d, bf_a_im_d;
    reg signed [WL-1:0] in_re_d, in_im_d;
    reg signed [WL-1:0] delay_out_re_d, delay_out_im_d;
    reg sel_d;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            bf_a_re_d <= 0; bf_a_im_d <= 0;
            in_re_d <= 0;   in_im_d <= 0;
            delay_out_re_d <= 0; delay_out_im_d <= 0;
            sel_d <= 0;
        end else if (valid_in) begin
            bf_a_re_d <= bf_a_re;
            bf_a_im_d <= bf_a_im;
            in_re_d   <= in_real;
            in_im_d   <= in_imag;
            delay_out_re_d <= delay_out_re;
            delay_out_im_d <= delay_out_im;
            sel_d <= sel;
        end
    end

    assign delay_in_re = sel_d ? mult_out_re : in_re_d;
    assign delay_in_im = sel_d ? mult_out_im : in_im_d;

    assign out_real = sel_d ? bf_a_re_d : delay_out_re_d;
    assign out_imag = sel_d ? bf_a_im_d : delay_out_im_d;
endmodule