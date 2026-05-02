module real_comparator #(

    parameter INPUT_WIDTH = 16,

    parameter NUM_SYMBOLS = 4096,

    parameter SYMBOL = 12

)(

    input wire clk,
    
    input wire rst_n,

    input wire valid_in,

    input wire signed [INPUT_WIDTH-1:0] fft_real_in,

    input wire signed [INPUT_WIDTH-1:0] fft_real_ref_in,

    input wire [SYMBOL-1:0] symbol_index,

    output reg comparison_result

);

    localparam THRESHOLD = NUM_SYMBOLS-20;   // majority

    reg [NUM_SYMBOLS-1:0] symbol_array;   
    //reg signed [INPUT_WIDTH-1:0] fft_real_in_reg;

    integer i;
    reg [$clog2(NUM_SYMBOLS):0] good_count;

    always @(posedge clk or negedge rst_n) begin

        if (!rst_n) begin
            // fft_real_in_reg <= 0;
            symbol_array <= {NUM_SYMBOLS{1'b0}};
            good_count <= 0;
            comparison_result <= 0;
        end
        else if (valid_in) begin
            // Register input
            // fft_real_in_reg <= fft_real_in;

            // Store comparison result
            symbol_array[symbol_index] <= 
                (fft_real_in[INPUT_WIDTH-1:0] == fft_real_ref_in[INPUT_WIDTH-1:0]);
            
            if(fft_real_in[INPUT_WIDTH-1:0] == fft_real_ref_in[INPUT_WIDTH-1:0]) begin
                good_count <= good_count + 1;
            end

            // Majority decision
            comparison_result <= (good_count >= THRESHOLD);
        end

    end

endmodule

//module real_comparator #(
//    parameter INPUT_WIDTH = 16,
//    parameter NUM_SYMBOLS = 4096,
//    parameter SYMBOL = 12
//)(
//    input wire clk,
//    input wire rst_n,
//    input wire valid_in,
//    input wire signed [INPUT_WIDTH-1:0] fft_real_in,
//    input wire signed [INPUT_WIDTH-1:0] fft_real_ref_in,
//    input wire [SYMBOL-1:0] symbol_index,
//    output reg comparison_result
//);

//    localparam THRESHOLD = NUM_SYMBOLS-20;   // majority

//    reg [NUM_SYMBOLS-1:0] symbol_array;   
//    reg signed [INPUT_WIDTH-1:0] fft_real_in_reg;
    
//    // NEW: Pipeline register to match the 1-cycle delay of fft_real_in_reg
//    reg valid_in_reg; 

//    reg [$clog2(NUM_SYMBOLS):0] good_count;

//    always @(posedge clk or negedge rst_n) begin
//        if (!rst_n) begin
//            fft_real_in_reg <= 0;
//            valid_in_reg <= 0;
//            symbol_array <= {NUM_SYMBOLS{1'b0}};
//            good_count <= 0;
//            comparison_result <= 0;
//        end
//        else begin
//            // 1. Pipeline the inputs
//            fft_real_in_reg <= fft_real_in;
//            valid_in_reg <= valid_in; // Now valid_in_reg aligns perfectly with fft_real_in_reg

//            // 2. ONLY evaluate and store when the aligned valid signal is high
//            if (valid_in_reg == 1'b1) begin
                
//                // Store comparison result
//                symbol_array[symbol_index] <= 
//                    (fft_real_in_reg[INPUT_WIDTH-1:0] == fft_real_ref_in[INPUT_WIDTH-1:0]);
                
//                // Increment score if they match
//                if(fft_real_in_reg[INPUT_WIDTH-1:0] == fft_real_ref_in[INPUT_WIDTH-1:0]) begin
//                    good_count <= good_count + 1; // Note: Changed to <= for sequential logic
//                end
                
//            end

//            // 3. Majority decision (only goes high if we are currently processing valid data)
//            comparison_result <= ((good_count >= THRESHOLD) && valid_in_reg);
//        end
//    end

//endmodule