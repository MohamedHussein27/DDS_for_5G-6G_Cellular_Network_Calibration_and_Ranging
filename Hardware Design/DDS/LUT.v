module LUT #(
    parameter memory_width = 8 ,
    parameter memory_depth = 16 
) (
    input clk,
    input rst_n,
    input [memory_depth-1 : 0] address, // Mapped address from quadrant logic
    input negative_flag,                // Instructs module to negate the output
    output reg [memory_width-1 :0] amplitude
);

// Memory array declaration. 
// Note: User must generate "memfile.mem" with quarter-wave sine data.
reg [memory_width-1 :0] memory [0 : (2**memory_depth)-1];

initial begin
    // Initialize memory from external hex file during synthesis/simulation
    $readmemh ("memfile.mem", memory);
end

always @(posedge clk, negedge rst_n) begin
    if (!rst_n) begin
        amplitude <= 0;
    end
    else begin
        // Synchronous read from memory
        if (negative_flag) begin
            // 3rd and 4th quadrants: Negate the amplitude.
            // Using 2's complement negation (-memory[address]). 
            // Ensure the downstream logic interprets this as a signed value.
            amplitude <= -memory[address];
        end
        else begin
            // 1st and 2nd quadrants: Output the positive amplitude directly.
            amplitude <= memory[address];
        end
    end
end
    
endmodule