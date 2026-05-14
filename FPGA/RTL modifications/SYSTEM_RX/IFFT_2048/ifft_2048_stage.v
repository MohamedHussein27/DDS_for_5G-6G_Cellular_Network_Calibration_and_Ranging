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


// module ifft_stage_2048 #(
//     parameter WL = 16,
//     parameter DELAY_LEN = 1024, 
//     parameter ROM_DEPTH = 1024
// )(
//     input wire clk,
//     input wire rst_n,
//     input wire valid_in, 
//     input wire sel,
//     input wire [(ROM_DEPTH == 1) ? 0 : $clog2(ROM_DEPTH)-1 : 0] addr,
//     input wire signed [WL-1:0] in_real,
//     input wire signed [WL-1:0] in_imag,
//     output wire signed [WL-1:0] out_real,
//     output wire signed [WL-1:0] out_imag
// );

//     // --- Internal Wire Declarations ---
//     wire signed [WL-1:0] delay_in_re, delay_in_im;
//     wire signed [WL-1:0] delay_out_re, delay_out_im;
//     wire signed [WL-1:0] bf_a_re, bf_a_im, bf_b_re, bf_b_im;
//     wire signed [WL-1:0] mult_out_re, mult_out_im;
    
//     // Wires for the raw ROM output
//     wire signed [WL-1:0] twiddle_re_rom, twiddle_im_rom;
//     // Wires for the conjugated twiddle used in multiplication
//     wire signed [WL-1:0] twiddle_re, twiddle_im;

//     // --- MUX Logic ---
//     assign delay_in_re = sel ? mult_out_re : in_real;
//     assign delay_in_im = sel ? mult_out_im : in_imag;
    
//     assign out_real = sel ? bf_a_re : delay_out_re;
//     assign out_imag = sel ? bf_a_im : delay_out_im;

//     // --- Sub-module Instantiations ---
    
//     // 1. Delay Feedback Buffer
//     delayfeedback #(
//         .WL(WL), 
//         .L(DELAY_LEN) 
//     ) delay_inst (
//         .clk(clk),
//         .rst_n(rst_n),
//         .en(valid_in),
//         .data_in_real(delay_in_re),
//         .data_in_imag(delay_in_im),
//         .data_out_real(delay_out_re),
//         .data_out_imag(delay_out_im)
//     );

//     // 2. Butterfly Unit
//     butterfly #(.WL(WL)) bf_inst (
//         .in1_real(delay_out_re), 
//         .in1_imag(delay_out_im),
//         .in2_real(in_real),      
//         .in2_imag(in_imag),
//         .a_real(bf_a_re),        
//         .a_imag(bf_a_im),
//         .b_real(bf_b_re),        
//         .b_imag(bf_b_im)
//     );

//     // 3. Twiddle Factor ROM
//     // Connect the 'addr' input and wire up the raw ROM outputs
//     twiddlerom_2048 #(
//         .WL(WL), 
//         .DEPTH(ROM_DEPTH)
//     ) rom_inst (
//         .clk(clk),
//         .addr(addr), // Fixed: Use the input addr instead of hardcoded 1'b0
//         .W_real(twiddle_re_rom), 
//         .W_img(twiddle_im_rom)
//     );

//     // 4. IFFT CONJUGATION
//     // To perform IFFT using an FFT architecture, we multiply by the conjugate 
//     // of the twiddle factor: W* = real - j(imag)
//     assign twiddle_re = twiddle_re_rom;
//     assign twiddle_im = -twiddle_im_rom; // Negate imaginary part for conjugation

//     // 5. Complex Multiplier
//     multiplier #(.WL(WL)) mult_inst (
//         .re1(bf_b_re),
//         .im1(bf_b_im),
//         .re2(twiddle_re), // Connected to conjugated twiddle
//         .im2(twiddle_im), // Connected to conjugated twiddle
//         .re_out(mult_out_re),
//         .im_out(mult_out_im)
//     );

// endmodule



module ifft_stage_2048 #(
    parameter WL = 16,
    parameter DELAY_LEN = 1024,
    parameter ROM_DEPTH = 1024
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in,
    input wire sel,

    input wire signed [WL-1:0] in_real,
    input wire signed [WL-1:0] in_imag,

    input wire signed [WL-1:0] twiddle_re,
    input wire signed [WL-1:0] twiddle_im,

    output wire signed [WL-1:0] out_real,
    output wire signed [WL-1:0] out_imag
);

    // =====================================================
    // Internal Signals
    // =====================================================
    wire signed [WL-1:0] delay_in_re, delay_in_im;
    wire signed [WL-1:0] delay_out_re, delay_out_im;

    wire signed [WL-1:0] bf_a_re, bf_a_im;
    wire signed [WL-1:0] bf_b_re, bf_b_im;

    wire signed [WL-1:0] mult_out_re, mult_out_im;

    wire signed [WL-1:0] twiddle_re_ifft, twiddle_im_ifft;

    // =====================================================
    // Delay Feedback
    // =====================================================
    delayfeedback #(
        .WL(WL),
        .L(DELAY_LEN - 1)
    ) delay_inst (
        .clk(clk),
        .rst_n(rst_n),
        .en(valid_in),

        .data_in_real(delay_in_re),
        .data_in_imag(delay_in_im),

        .data_out_real(delay_out_re),
        .data_out_imag(delay_out_im)
    );

    // =====================================================
    // Butterfly
    // =====================================================
    butterfly #(.WL(WL)) bf_inst (
        .in1_real(delay_out_re),
        .in1_imag(delay_out_im),

        .in2_real(in_real),
        .in2_imag(in_imag),

        .a_real(bf_a_re),
        .a_imag(bf_a_im),

        .b_real(bf_b_re),
        .b_imag(bf_b_im)
    );

    // =====================================================
    // Twiddle ROM
    // =====================================================
    // twiddlerom_2048 #(
    //     .WL(WL),
    //     .DEPTH(ROM_DEPTH)
    // ) rom_inst (
    //     .addr(addr),
    //     .clk(clk),

    //     .W_real(twiddle_re),
    //     .W_img(twiddle_im)
    // );

    // =====================================================
    // IFFT Conjugation
    // =====================================================
    assign twiddle_re_ifft = twiddle_re;
    assign twiddle_im_ifft = -twiddle_im;

    // =====================================================
    // Complex Multiplier
    // =====================================================
    multiplier #(.WL(WL)) mult_inst (
        .clk(clk),
        .rst_n(rst_n),
        .en(valid_in),

        .re1(bf_b_re),
        .im1(bf_b_im),

        .re2(twiddle_re_ifft),
        .im2(twiddle_im_ifft),

        .re_out(mult_out_re),
        .im_out(mult_out_im)
    );

    // =====================================================
    // Pipeline Alignment Registers
    // =====================================================
    reg signed [WL-1:0] bf_a_re_d, bf_a_im_d;
    reg signed [WL-1:0] in_re_d, in_im_d;
    reg signed [WL-1:0] delay_out_re_d, delay_out_im_d;

    reg sel_d;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            bf_a_re_d <= 0;
            bf_a_im_d <= 0;

            in_re_d <= 0;
            in_im_d <= 0;

            delay_out_re_d <= 0;
            delay_out_im_d <= 0;

            sel_d <= 0;
        end
        else if (valid_in) begin
            bf_a_re_d <= bf_a_re;
            bf_a_im_d <= bf_a_im;

            in_re_d <= in_real;
            in_im_d <= in_imag;

            delay_out_re_d <= delay_out_re;
            delay_out_im_d <= delay_out_im;

            sel_d <= sel;
        end
    end

    // =====================================================
    // Final Routing
    // =====================================================
    assign delay_in_re = sel_d ? mult_out_re : in_re_d;
    assign delay_in_im = sel_d ? mult_out_im : in_im_d;

    assign out_real = sel_d ? bf_a_re_d : delay_out_re_d;
    assign out_imag = sel_d ? bf_a_im_d : delay_out_im_d;

endmodule