// Analog Devices 
// GP Ain-shams University
// Butterfly unit for 4096-FFT 


/* Description ..........
    This module now drives the pipeline by counting every clock cycle.
*/

module controlunit_4096 #(
    parameter N = 4096 
) (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        valid_in,
    output wire [11:0] addr, 
    output wire [11:0] sel,
    output wire        pipeline_en, 
    output wire        valid_out
);

    // --------------------------------------------------
    // Internal Signals
    // --------------------------------------------------
    reg [11:0] count;
    reg [N-2:0] valid_pipe;

    // --------------------------------------------------
    // Pipeline Enable (active when data is flowing)
    // --------------------------------------------------
    assign pipeline_en = valid_in | valid_out;

    // --------------------------------------------------
    // Sequential Logic
    // --------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            count      <= 12'd0;
            valid_pipe <= { (N-1){1'b0} };
        end else begin
            // Shift valid signal through pipeline (N-1 delay)
            valid_pipe <= {valid_pipe[N-3:0], valid_in};

            // Count ONLY when pipeline is active
            if (pipeline_en)
                count <= count + 1'b1;
            else
                count <= 12'd0;   // reset when idle (important for frame alignment)
        end
    end

    // --------------------------------------------------
    // Outputs
    // --------------------------------------------------
    assign addr = count;
    assign sel  = count;

    // Output valid after N-1 cycles delay
    assign valid_out = valid_pipe[N-2];

endmodule