module simple_dut (
    input  logic        clk,
    input  logic        rst_n,
    input  logic        valid,
    input  logic [7:0]  data_in,
    output logic [7:0]  data_out,
    output logic        ready
);
//(multiply input by 2 with 2-cycle delay).
    logic [7:0] stage1, stage2;
    logic       valid_d1, valid_d2;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            stage1   <= 0;
            stage2   <= 0;
            valid_d1 <= 0;
            valid_d2 <= 0;
        end else begin
            stage1   <= data_in;
            valid_d1 <= valid;

            stage2   <= stage1 * 2;
            valid_d2 <= valid_d1;
        end
    end

    assign data_out = stage2;
    assign ready    = valid_d2;

endmodule