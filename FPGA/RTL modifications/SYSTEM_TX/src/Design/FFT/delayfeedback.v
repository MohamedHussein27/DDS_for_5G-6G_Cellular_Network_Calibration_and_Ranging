module delayfeedback #(
    parameter L = 2048,
    parameter WL = 16
)(
    input wire clk,
    input wire rst_n,
    input wire en,
    input wire signed [WL-1:0] data_in_real,
    input wire signed [WL-1:0] data_in_imag,
    output wire signed [WL-1:0] data_out_real,
    output wire signed [WL-1:0] data_out_imag
);
    generate
        if (L == 0) begin : gen_wire
            assign data_out_real = data_in_real;
            assign data_out_imag = data_in_imag;
        end else begin : gen_shift
            localparam ARR_MAX = (L > 0) ? (L - 1) : 0;
            reg signed [WL-1:0] delay_reg_real [0:ARR_MAX];
            reg signed [WL-1:0] delay_reg_imag [0:ARR_MAX];
            integer i, j;

            always @(posedge clk or negedge rst_n) begin 
                if (!rst_n) begin
                    for (i = 0; i <= ARR_MAX; i = i + 1) begin
                        delay_reg_real[i] <= 0;
                        delay_reg_imag[i] <= 0;
                    end
                end else if (en) begin
                    for (j = ARR_MAX; j > 0; j = j - 1) begin
                        delay_reg_real[j] <= delay_reg_real[j-1];
                        delay_reg_imag[j] <= delay_reg_imag[j-1];
                    end
                    delay_reg_real[0] <= data_in_real;
                    delay_reg_imag[0] <= data_in_imag;
                end
            end

            assign data_out_real = delay_reg_real[ARR_MAX];
            assign data_out_imag = delay_reg_imag[ARR_MAX];
        end
    endgenerate
endmodule

