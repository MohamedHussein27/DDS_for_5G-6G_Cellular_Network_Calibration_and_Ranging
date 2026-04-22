module first_quad_address #(
    parameter address_width = 16,
    parameter quartre_memory = 0.25*(2**(address_width)),
    parameter half_memory = 0.5*(2**(address_width)),
    parameter three_quartre_memory = 0.75*(2**(address_width))
) (
    input [address_width-1:0] truncated_address,
    input clk,
    output reg [address_width-1:0] first_address,
    output reg negative_flag
);



always @(*) begin

        case (truncated_address[address_width-1:address_width-2])
        2'b11 : begin         //4th quad
            first_address = (2**address_width)-1-truncated_address;
          //  negative_flag = 1'b1;
        end
        2'b10 : begin        //3rd quad
            first_address = truncated_address - ( half_memory );
           // negative_flag = 1'b1;
        end
        2'b01 : begin        //2nd quad
            first_address = half_memory - 1 - truncated_address;
        //    negative_flag = 1'b0;
        end
        2'b00 : begin        //1st quad
            first_address = truncated_address;
           // negative_flag = 1'b0;
        end
        endcase
end


always @(posedge clk) begin

        case (truncated_address[address_width-1:address_width-2])
        2'b11 : begin         //4th quad
           // first_address = (2**address_width)-1-truncated_address;
            negative_flag <= 1'b1;
        end
        2'b10 : begin        //3rd quad
         //   first_address <= truncated_address - ( half_memory );
            negative_flag <= 1'b1;
        end
        2'b01 : begin        //2nd quad
         //   first_address <= half_memory - 1 - truncated_address;
            negative_flag <= 1'b0;
        end
        2'b00 : begin        //1st quad
         //   first_address <= truncated_address;
            negative_flag <= 1'b0;
        end
        endcase
end


endmodule