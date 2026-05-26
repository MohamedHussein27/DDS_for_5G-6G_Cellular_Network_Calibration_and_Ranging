module imaginary_comparator #(

    parameter INPUT_WIDTH = 16,

    parameter NUM_SYMBOLS = 4096,

    parameter SYMBOL = 12                                 //???????????????????

)(

    input wire clk,
    
    input wire rst_n,

    input wire ofdm_valid_in,

    input wire radar_valid_in,

    input wire signed [INPUT_WIDTH-1:0] OFDM_imag_in,

    input wire signed [INPUT_WIDTH-1:0] OFDM_imag_ref_in,

    input wire signed [INPUT_WIDTH-1:0] radar_imag_in,

    input wire signed [INPUT_WIDTH-1:0] radar_imag_ref_in,

    input wire [SYMBOL-1:0] symbol_index_ofdm,

    input wire [SYMBOL-1:0] symbol_index_radar,

    output reg comparison_result_ofdm,

    output reg comparison_result_radar

);

    localparam THRESHOLD = NUM_SYMBOLS-500;   // majority

    // reg [NUM_SYMBOLS-1:0] symbol_array_ofdm;   
    reg signed [INPUT_WIDTH-1:0] ofdm_imag_in_reg;

    integer i;
    reg [$clog2(NUM_SYMBOLS):0] good_count_ofdm;

    always @(posedge clk or negedge rst_n) begin

        if (!rst_n) begin
            ofdm_imag_in_reg <= 0;
            // symbol_array_ofdm <= {NUM_SYMBOLS{1'b0}};
            good_count_ofdm <= 0;
            comparison_result_ofdm <= 0;
        end
        else if (ofdm_valid_in) begin
            // Register input
            ofdm_imag_in_reg <= OFDM_imag_in;

            // // Store comparison result
            // symbol_array_ofdm[symbol_index_ofdm] <= 
            //     (ofdm_imag_in_reg[INPUT_WIDTH-1:6] == OFDM_imag_ref_in[INPUT_WIDTH-1:6]);
            
            if(ofdm_imag_in_reg[INPUT_WIDTH-1:6] == OFDM_imag_ref_in[INPUT_WIDTH-1:6]) begin
                good_count_ofdm <= good_count_ofdm + 1;
            end

            // Majority decision
            comparison_result_ofdm <= (good_count_ofdm >= THRESHOLD);
        end

    end

    // reg [NUM_SYMBOLS-1:0] symbol_array_radar;   
    reg signed [INPUT_WIDTH-1:0] radar_imag_in_reg;

    reg [$clog2(NUM_SYMBOLS):0] good_count_radar;

    always @(posedge clk or negedge rst_n) begin

        if (!rst_n) begin
            radar_imag_in_reg <= 0;
            // symbol_array_radar <= {NUM_SYMBOLS{1'b0}};
            good_count_radar <= 0;
            comparison_result_radar <= 0;
        end
        else if (radar_valid_in) begin
            // Register input
            radar_imag_in_reg <= radar_imag_in;

            // Store comparison result
            // symbol_array_radar[symbol_index_radar] <= 
            //     (radar_imag_in_reg[INPUT_WIDTH-1:6] == radar_imag_ref_in[INPUT_WIDTH-1:6]);

            if(radar_imag_in_reg[INPUT_WIDTH-1:6] == radar_imag_ref_in[INPUT_WIDTH-1:6]) begin
                good_count_radar <= good_count_radar + 1;
            end

            // Majority decision
            comparison_result_radar <= (good_count_radar >= THRESHOLD);
        end

    end

endmodule
