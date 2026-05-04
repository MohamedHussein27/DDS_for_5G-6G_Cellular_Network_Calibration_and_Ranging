//module comparator #(

//    parameter INPUT_WIDTH = 8,

//    parameter NUM_SYMBOLS = 4096,

//    parameter SYMBOL = 12

//)(

//    input wire clk,
    
//    input wire rst_n,

//    input wire signed [INPUT_WIDTH-1:0] dds_in,

//    input wire signed [INPUT_WIDTH-1:0] ref_in,

//    input wire [SYMBOL-1:0] symbol_index,

//    output reg comparison_result

//);

//    localparam THRESHOLD = NUM_SYMBOLS/2;   // majority

//    reg [NUM_SYMBOLS-1:0] symbol_array;   
//    reg signed [INPUT_WIDTH-1:0] dds_in_reg;

//    integer i;
//    reg [$clog2(NUM_SYMBOLS):0] good_count;

//    always @(posedge clk or negedge rst_n) begin
//        if (!rst_n) begin
//            dds_in_reg <= 0;
//            symbol_array <= {NUM_SYMBOLS{1'b0}};
//            good_count <= 0;
//            comparison_result <= 0;
//        end
//        else begin

//            // Register input
//            dds_in_reg <= dds_in;

//            // Store comparison result
//            symbol_array[symbol_index] <= 
//                (dds_in_reg[INPUT_WIDTH-1:0] == ref_in[INPUT_WIDTH-1:0]);

//            // Count number of '1's
//            good_count = 0;
//            for (i = 0; i < NUM_SYMBOLS; i = i + 1) begin
//                good_count = good_count + symbol_array[i];
//            end

//            // Majority decision
//            comparison_result <= (good_count >= THRESHOLD);

//        end
//    end

//endmodule

module imaginary_comparator #(

    parameter INPUT_WIDTH = 16,

    parameter NUM_SYMBOLS = 4096,

    parameter SYMBOL = 12                                 //???????????????????

)(

    input wire clk,
    
    input wire rst_n,

    input wire valid_in,

    input wire signed [INPUT_WIDTH-1:0] fft_imag_in,

    input wire signed [INPUT_WIDTH-1:0] fft_imag_ref_in,

    input wire [SYMBOL-1:0] symbol_index,

    output reg comparison_result

);

    localparam THRESHOLD = NUM_SYMBOLS-500;   // majority

    reg [NUM_SYMBOLS-1:0] symbol_array;   
    reg signed [INPUT_WIDTH-1:0] fft_imag_in_reg;

    integer i;
    reg [$clog2(NUM_SYMBOLS):0] good_count;

    always @(posedge clk or negedge rst_n) begin

        if (!rst_n) begin
            fft_imag_in_reg <= 0;
            symbol_array <= {NUM_SYMBOLS{1'b0}};
            good_count <= 0;
            comparison_result <= 0;
        end
        else if (valid_in) begin
            // Register input
            fft_imag_in_reg <= fft_imag_in;

            // Store comparison result
            symbol_array[symbol_index] <= 
                (fft_imag_in_reg[INPUT_WIDTH-1:6] == fft_imag_ref_in[INPUT_WIDTH-1:6]);
            
            if(fft_imag_in_reg[INPUT_WIDTH-1:6] == fft_imag_ref_in[INPUT_WIDTH-1:6]) begin
                good_count <= good_count + 1;
            end

            // Majority decision
            comparison_result <= (good_count >= THRESHOLD);
        end

    end

endmodule