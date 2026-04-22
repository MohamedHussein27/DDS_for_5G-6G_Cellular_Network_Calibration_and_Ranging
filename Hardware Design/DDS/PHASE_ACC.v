module PHASE_ACC #(
    parameter TUNING_WORD_WIDTH = 32,
    parameter CYCLES_WIDTH = 13,
    parameter F_START_WIDTH = 3,
    parameter TRUNCATED_ADDRESS_WIDTH = 16  
) (
    input clk,
    input rst_n,
    input [TUNING_WORD_WIDTH-1:0] FTW_start,
    input [CYCLES_WIDTH-1:0] cycles,
    input [TUNING_WORD_WIDTH-1:0] FTW_step,
    output [TRUNCATED_ADDRESS_WIDTH-1:0] truncated_address
);

reg [CYCLES_WIDTH-1:0] cycles_counter;
reg [TUNING_WORD_WIDTH-1:0] current_ftw;
reg [TUNING_WORD_WIDTH-1:0] truncated_phase;

always @(posedge clk, negedge rst_n) begin
    if (!rst_n) begin
        cycles_counter <=0;
        current_ftw <=FTW_start;
        truncated_phase <=0;
    end

    else begin

        if (cycles_counter==cycles+1) begin
        cycles_counter <=0;
        current_ftw <=FTW_start;
        truncated_phase <=0;
        end

        else begin
            current_ftw <= current_ftw + FTW_step;
            truncated_phase <= truncated_phase + current_ftw;
            cycles_counter <= cycles_counter + 1;

        end

    end
end

assign truncated_address = truncated_phase [TUNING_WORD_WIDTH-1 : TUNING_WORD_WIDTH-1-(TRUNCATED_ADDRESS_WIDTH-1) ] ;





endmodule