module bram_write_counter #(
    parameter ADDR_WIDTH = 11
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in_ofdm,
    input wire valid_in_radar,
    output reg [ADDR_WIDTH-1:0] write_addr_ofdm,
    output reg [ADDR_WIDTH-1:0] write_addr_radar
);

    always @(posedge clk) begin
        if (!rst_n) begin
            write_addr_ofdm <= 0;
        end else if (valid_in_ofdm)
                                 begin
            // Increment address only when DDS data is valid
            // Stop at the max value to prevent overwriting (or remove the if statement to let it loop)
            if (write_addr_ofdm < (2**ADDR_WIDTH) - 1) begin
                write_addr_ofdm <= write_addr_ofdm + 1;
            end
        end
    end
    
    always @(posedge clk) begin
        if (!rst_n) begin
            write_addr_radar <= 0;
        end else if (valid_in_radar)
                                 begin
            // Increment address only when DDS data is valid
            // Stop at the max value to prevent overwriting (or remove the if statement to let it loop)
            if (write_addr_radar < (2**ADDR_WIDTH) - 1) begin
                write_addr_radar <= write_addr_radar + 1;
            end
        end
    end
endmodule