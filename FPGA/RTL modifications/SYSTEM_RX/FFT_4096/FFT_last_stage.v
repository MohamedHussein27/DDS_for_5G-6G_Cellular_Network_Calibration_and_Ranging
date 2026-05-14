module fft_stage_last #(
    parameter WL = 16
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in, 
    input wire sel,
    input wire signed [WL-1:0] in_real,
    input wire signed [WL-1:0] in_imag,
    output reg signed [WL-1:0] out_real,
    output reg signed [WL-1:0] out_imag
);
    // Stage 12 replaces the Multiplier and ROM with direct Butterfly routing.
    // It maintains the exact 1-clock-cycle latency required by the pipeline.
    reg signed [WL-1:0] delay_reg_re, delay_reg_im;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            out_real     <= 0; 
            out_imag     <= 0;
            delay_reg_re <= 0; 
            delay_reg_im <= 0;
        end else if (valid_in) begin
            if (!sel) begin
                // Feed-through phase
                out_real     <= delay_reg_re;
                out_imag     <= delay_reg_im;
                delay_reg_re <= in_real;
                delay_reg_im <= in_imag;
            end else begin
                // Butterfly phase (W = 1, so no multiplication needed)
                out_real     <= delay_reg_re + in_real;
                out_imag     <= delay_reg_im + in_imag;
                delay_reg_re <= delay_reg_re - in_real;
                delay_reg_im <= delay_reg_im - in_imag;
            end
        end
    end
endmodule