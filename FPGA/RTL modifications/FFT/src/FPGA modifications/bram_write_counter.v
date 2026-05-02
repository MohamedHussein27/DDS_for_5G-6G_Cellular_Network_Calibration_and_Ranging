module bram_write_counter #(
    parameter ADDR_WIDTH = 12
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in,
    output reg [ADDR_WIDTH-1:0] write_addr
);
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            write_addr <= 0;
        end else if (valid_in)
                                 begin
            // Increment address only when DDS data is valid
            // Stop at the max value to prevent overwriting (or remove the if statement to let it loop)
            if (write_addr < (2**ADDR_WIDTH) - 1) begin
                write_addr <= write_addr + 1;
            end
        end
    end
    
endmodule