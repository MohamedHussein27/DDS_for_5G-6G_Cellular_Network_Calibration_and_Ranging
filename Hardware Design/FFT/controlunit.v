// Analog Devices 
// GP Ain-shams University
// Butterfly unit for 4096-FFT 


/* Description ..........
    This module now drives the pipeline by counting every clock cycle.
*/

module controlunit #(
    parameter N = 4096 // FFT size
) (
    input wire clk,
    input wire rst_n,
    output wire [11:0] addr, 
    output wire [11:0] sel 
);
    
    reg [11:0] count;

    // 12-bit Counter
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            count <= 12'd0;
        end else begin
            count <= count + 1'b1;
        end
    end

    // The counter bits directly drive the selects and addresses
    assign sel = count;
    assign addr = count;

endmodule