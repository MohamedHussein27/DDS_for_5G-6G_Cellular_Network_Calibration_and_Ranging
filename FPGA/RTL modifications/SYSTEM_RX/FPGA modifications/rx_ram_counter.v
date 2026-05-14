module rx_ram_counter #(
    parameter ADDR_WIDTH = 11
)(
    input wire clk,
    input wire rst_n,
    input wire en,
    output reg [ADDR_WIDTH-1:0] count
);
    always @(posedge clk) begin
        if (!rst_n) begin
            count <= 0;
        end else if (en) begin
            count <= count + 1;
        end
    end
endmodule
