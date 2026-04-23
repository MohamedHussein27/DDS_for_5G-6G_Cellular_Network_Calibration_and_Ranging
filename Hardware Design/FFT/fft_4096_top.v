module fft_4096_top #(
    parameter WL = 16
)(
    input wire clk,
    input wire rst_n,
    input wire signed [WL-1:0] in_real,
    input wire signed [WL-1:0] in_imag,
    output wire signed [WL-1:0] out_real,
    output wire signed [WL-1:0] out_imag
);

    // Global control signals
    wire [11:0] global_addr;
    wire [11:0] global_sel;

    controlunit #(.N(4096)) ctrl (
        .clk(clk),
        .rst_n(rst_n),
        .addr(global_addr),
        .sel(global_sel)
    );

    // Interconnect wires between stages
    // node[0] is the top-level input, node[12] is the top-level output
    wire signed [WL-1:0] stage_re [0:12];
    wire signed [WL-1:0] stage_im [0:12];

    assign stage_re[0] = in_real;
    assign stage_im[0] = in_imag;

    assign out_real = stage_re[12];
    assign out_imag = stage_im[12];

    // Stage 1 (Delay = 2048)
    fft_stage #(.WL(WL), .DELAY_LEN(2048), .ROM_DEPTH(2048)) stg1 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[11]), .addr(global_addr[10:0]),
        .in_real(stage_re[0]), .in_imag(stage_im[0]), .out_real(stage_re[1]), .out_imag(stage_im[1])
    );

    // Stage 2 (Delay = 1024)
    fft_stage #(.WL(WL), .DELAY_LEN(1024), .ROM_DEPTH(1024)) stg2 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[10]), .addr(global_addr[9:0]),
        .in_real(stage_re[1]), .in_imag(stage_im[1]), .out_real(stage_re[2]), .out_imag(stage_im[2])
    );

    // Stage 3 (Delay = 512)
    fft_stage #(.WL(WL), .DELAY_LEN(512), .ROM_DEPTH(512)) stg3 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[9]), .addr(global_addr[8:0]),
        .in_real(stage_re[2]), .in_imag(stage_im[2]), .out_real(stage_re[3]), .out_imag(stage_im[3])
    );

    // Stage 4 (Delay = 256)
    fft_stage #(.WL(WL), .DELAY_LEN(256), .ROM_DEPTH(256)) stg4 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[8]), .addr(global_addr[7:0]),
        .in_real(stage_re[3]), .in_imag(stage_im[3]), .out_real(stage_re[4]), .out_imag(stage_im[4])
    );

    // Stage 5 (Delay = 128)
    fft_stage #(.WL(WL), .DELAY_LEN(128), .ROM_DEPTH(128)) stg5 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[7]), .addr(global_addr[6:0]),
        .in_real(stage_re[4]), .in_imag(stage_im[4]), .out_real(stage_re[5]), .out_imag(stage_im[5])
    );

    // Stage 6 (Delay = 64)
    fft_stage #(.WL(WL), .DELAY_LEN(64), .ROM_DEPTH(64)) stg6 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[6]), .addr(global_addr[5:0]),
        .in_real(stage_re[5]), .in_imag(stage_im[5]), .out_real(stage_re[6]), .out_imag(stage_im[6])
    );

    // Stage 7 (Delay = 32)
    fft_stage #(.WL(WL), .DELAY_LEN(32), .ROM_DEPTH(32)) stg7 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[5]), .addr(global_addr[4:0]),
        .in_real(stage_re[6]), .in_imag(stage_im[6]), .out_real(stage_re[7]), .out_imag(stage_im[7])
    );

    // Stage 8 (Delay = 16)
    fft_stage #(.WL(WL), .DELAY_LEN(16), .ROM_DEPTH(16)) stg8 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[4]), .addr(global_addr[3:0]),
        .in_real(stage_re[7]), .in_imag(stage_im[7]), .out_real(stage_re[8]), .out_imag(stage_im[8])
    );

    // Stage 9 (Delay = 8)
    fft_stage #(.WL(WL), .DELAY_LEN(8), .ROM_DEPTH(8)) stg9 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[3]), .addr(global_addr[2:0]),
        .in_real(stage_re[8]), .in_imag(stage_im[8]), .out_real(stage_re[9]), .out_imag(stage_im[9])
    );

    // Stage 10 (Delay = 4)
    fft_stage #(.WL(WL), .DELAY_LEN(4), .ROM_DEPTH(4)) stg10 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[2]), .addr(global_addr[1:0]),
        .in_real(stage_re[9]), .in_imag(stage_im[9]), .out_real(stage_re[10]), .out_imag(stage_im[10])
    );

    // Stage 11 (Delay = 2)
    fft_stage #(.WL(WL), .DELAY_LEN(2), .ROM_DEPTH(2)) stg11 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[1]), .addr(global_addr[0]),
        .in_real(stage_re[10]), .in_imag(stage_im[10]), .out_real(stage_re[11]), .out_imag(stage_im[11])
    );

    // Stage 12 (Delay = 1)
    // Note: The twiddle factor for a 2-point FFT is always 1 (W_2^0).
    fft_stage #(.WL(WL), .DELAY_LEN(1), .ROM_DEPTH(1)) stg12 (
        .clk(clk), .rst_n(rst_n),
        .sel(global_sel[0]), .addr(1'b0), // Address is tied to 0 since ROM depth is 1
        .in_real(stage_re[11]), .in_imag(stage_im[11]), .out_real(stage_re[12]), .out_imag(stage_im[12])
    );

endmodule