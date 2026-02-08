// ---------------------------------------------------------
// Design: Simple Sequence Detector (Detects "101")
// ---------------------------------------------------------
module seq_detector (
    input  wire clk,
    input  wire rst_n,
    input  wire data_in,
    output reg  seq_detected
);

    // States
    localparam IDLE = 2'b00;
    localparam S1   = 2'b01; 
    localparam S10  = 2'b10; 
    localparam S101 = 2'b11; 

    reg [1:0] current_state, next_state;

    // State Memory
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) current_state <= IDLE;
        else        current_state <= next_state;
    end

    // Next State Logic
    always @(*) begin
        case (current_state)
            IDLE: next_state = (data_in) ? S1 : IDLE;
            S1:   next_state = (data_in) ? S1 : S10;
            S10:  next_state = (data_in) ? S101 : IDLE;
            S101: next_state = (data_in) ? S1 : IDLE; // Overlapping detection
            default: next_state = IDLE;
        endcase
    end

    // Output Logic
    always @(*) begin
        seq_detected = (current_state == S101);
    end

endmodule