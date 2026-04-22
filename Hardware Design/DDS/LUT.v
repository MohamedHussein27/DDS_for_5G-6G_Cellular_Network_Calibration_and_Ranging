module LUT #(
    parameter memory_width = 8 ,
    parameter memory_depth = 16 
) (
    input clk,
    input rst_n,
    input [memory_depth-1 : 0] address,
    
    output reg [memory_width-1 :0] amplitude
);

reg [memory_width-1 :0] memory [0 : (2**memory_depth)-1] ;

initial begin
    $readmemh ("memfile.mem", memory);
end

always @(posedge clk, negedge rst_n) begin
    if (!rst_n) begin
        amplitude <= 0;
    end

    else begin
            amplitude <= memory[address];
    end
end
    
endmodule