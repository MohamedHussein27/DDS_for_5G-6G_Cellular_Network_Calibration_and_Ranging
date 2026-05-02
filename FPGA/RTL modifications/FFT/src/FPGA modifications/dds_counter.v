module dds_counter #(
  //  parameter CYCLES = 4096,
    parameter ADDR_WIDTH = 12
) (
    input wire enable,
    input wire clk,
    input wire rst_n,
    output reg [ADDR_WIDTH-1:0] count
);
    always @(posedge clk , negedge rst_n) begin
        if (!rst_n) begin
            count <= 2;
        end

        else if (enable) begin
            count <= count + 1;
        end
    end
endmodule