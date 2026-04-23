// Analog Devices 
// GP Ain-shams University
// Butterfly unit for 4096-FFT 


/* Description ..........
    This module now drives the addr by counting every clock cycle.
    and controlling the valid out by delaying the valid in by N cycles.
*/

module controlunit #(
    parameter N = 4096 
) (
    input wire clk,
    input wire rst_n,
    input wire valid_in,     // High when DDS starts sending
    output wire [11:0] addr, 
    output wire [11:0] sel,
    output wire pipeline_en, 
    output wire valid_out    
);
    reg [11:0] count;
    // Shift register to delay the valid signal by N-1 cycles
    reg [N-2:0] valid_pipe; 

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            count <= 12'd0;
            valid_pipe <= 0;
        end else begin
            count <= count + 1'b1; 
            // Standard SDF delay: Output valid is input valid delayed by N-1
            valid_pipe <= {valid_pipe[N-3:0], valid_in}; 
        end
    end

    assign sel = count;
    assign addr = count;
    assign pipeline_en = 1'b1; // Always running to allow overlap shifting values whether it is input or output 
    assign valid_out = valid_pipe[N-2]; 
endmodule
