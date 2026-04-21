module first_quad_address #(
    parameter address_width = 16,
    parameter quartre_memory = 0.25*(2**(address_width)),
    parameter half_memory = 0.5*(2**(address_width)),
    parameter three_quartre_memory = 0.75*(2**(address_width))
) (
    input [address_width-1:0] truncated_address,
    output reg [address_width-1:0] first_address, // Address folded into [0, pi/2]
    output reg negative_flag                      // 1 = invert amplitude (pi to 2pi)
);

always @(*) begin
    // The top two bits divide the circle into 4 quadrants
    case (truncated_address[address_width-1:address_width-2])
        
        // 4th Quadrant (270 to 360 degrees)
        // Sine wave is negative and rising back to 0. We read the LUT backwards.
        2'b11 : begin         
            first_address = (2**address_width) - 1 - truncated_address;
            negative_flag = 1'b1;
        end
        
        // 3rd Quadrant (180 to 270 degrees)
        // Sine wave is negative and falling. We read the LUT forwards.
        2'b10 : begin        
            first_address = truncated_address - half_memory;
            negative_flag = 1'b1;
        end
        
        // 2nd Quadrant (90 to 180 degrees)
        // Sine wave is positive and falling. We read the LUT backwards.
        2'b01 : begin        
            first_address = half_memory - 1 - truncated_address;
            negative_flag = 1'b0;
        end
        
        // 1st Quadrant (0 to 90 degrees)
        // Sine wave is positive and rising. We read the LUT forwards directly.
        2'b00 : begin        
            first_address = truncated_address;
            negative_flag = 1'b0;
        end
        
    endcase
end

endmodule