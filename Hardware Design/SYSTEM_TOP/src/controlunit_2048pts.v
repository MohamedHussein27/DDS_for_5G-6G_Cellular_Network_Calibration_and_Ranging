// Analog Devices 
// GP Ain-shams University
// Butterfly unit for 2048-FFT 


/* Description ..........
    This module now drives the pipeline by counting every clock cycle.
*/

module controlunit_2048 #(
    parameter N = 2048 
) (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        valid_in,
    output wire [10:0] addr, 
    output wire [10:0] sel,
    output wire        pipeline_en, 
    output wire        valid_out
);

    // --------------------------------------------------
    // Internal Signals
    // --------------------------------------------------
    reg [10:0] count;
    reg [N-2:0] valid_pipe;

    // --------------------------------------------------
    // Output valid (delay = N-1 cycles)
    // --------------------------------------------------
    assign valid_out = valid_pipe[N-2];

    // --------------------------------------------------
    // Pipeline Enable (active when data flows)
    // --------------------------------------------------
    assign pipeline_en = valid_in | valid_out;

    // --------------------------------------------------
    // Sequential Logic
    // --------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            count      <= 11'd0;
            valid_pipe <= { (N-1){1'b0} };
        end else begin
            // Shift valid signal through pipeline
            valid_pipe <= {valid_pipe[N-3:0], valid_in};

            // Counter follows data movement
            if (pipeline_en)
                count <= count + 1'b1;
            else
                count <= 11'd0;   // reset when idle (frame alignment)
        end
    end

    // --------------------------------------------------
    // Outputs
    // --------------------------------------------------
    assign addr = count;
    assign sel  = count;

endmodule