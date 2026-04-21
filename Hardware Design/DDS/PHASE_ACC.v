module PHASE_ACC #(
    parameter TUNING_WORD_WIDTH = 32,
    parameter CYCLES_WIDTH = 13,
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
        // Asynchronous reset: Clear counters and reset to initial frequency
        cycles_counter <= 0;
        current_ftw    <= FTW_start;
        truncated_phase<= 0;
    end
    else begin
        // Check if the programmed duration has elapsed
        if (cycles_counter == cycles + 1) begin
            // Reset sequence for the next sweep/burst
            cycles_counter <= 0;
            current_ftw    <= FTW_start;
            truncated_phase<= 0;
        end
        else begin
            // Accumulate the Tuning Word to step the frequency (Chirp generation)
            current_ftw     <= current_ftw + FTW_step;
            
            // Accumulate the phase based on the current frequency
            // Overflows naturally, creating the sawtooth phase wrap-around
            truncated_phase <= truncated_phase + current_ftw;
            
            // Track how many cycles have passed
            cycles_counter  <= cycles_counter + 1;
        end
    end
end

// Phase Truncation: 
// Discard the lower bits of the accumulator to match the LUT address size.
// Keeping the MSBs preserves the structural phase information while discarding 
// fine resolution bits that don't fit in the LUT. This introduces phase truncation error (spurs).
assign truncated_address = truncated_phase [TUNING_WORD_WIDTH-1 : TUNING_WORD_WIDTH-1-(TRUNCATED_ADDRESS_WIDTH-1)];

endmodule