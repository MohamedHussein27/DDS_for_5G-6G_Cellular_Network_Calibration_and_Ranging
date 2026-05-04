module dds_counter #(
  //  parameter CYCLES = 4096,
    parameter ADDR_WIDTH = 2
) (
    input wire enable,
    input wire clk,
    input wire rst_n,
    output reg [ADDR_WIDTH-1:0] count
//    output reg enable_out
);
    always @(posedge clk , negedge rst_n) begin
        if (!rst_n) begin
            count <= 0;
//            enable_out <= 0;
        end else if (enable) begin
            if(count < 2**ADDR_WIDTH - 1) begin
                count <= count + 1;
//                enable_out <= 0;
//            end else begin
//                enable_out <= 1;
            end
        end
    end
endmodule