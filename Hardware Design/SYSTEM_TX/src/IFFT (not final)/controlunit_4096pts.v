// Analog Devices 
// GP Ain-shams University
// Butterfly unit for 4096-FFT 


/* Description ..........
    This module now drives the pipeline by counting every clock cycle.
*/

// module controlunit_4096 #(
//     parameter N = 4096 
// ) (
//     input wire clk,
//     input wire rst_n,
//     input wire valid_in,     // High when DDS starts sending
//     output wire [11:0] addr, 
//     output wire [11:0] sel,
//     output wire pipeline_en, 
//     output wire valid_out    
// );
//     reg [11:0] count;
//     // Shift register to delay the valid signal by N-1 cycles
//     reg [N-2:0] valid_pipe; 

//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             count <= 12'd0;
//             valid_pipe <= 0;
//         end else begin
//             count <= count + 1'b1; 
//             // Standard SDF delay: Output valid is input valid delayed by N-1
//             valid_pipe <= {valid_pipe[N-3:0], valid_in}; 
//         end
//     end

//     assign sel = count;
//     assign addr = count;
//     assign pipeline_en = 1'b1; // Always running to allow overlap
//     assign valid_out = valid_pipe[N-2]; 
// endmodule

// reg [11:0] count;        // enough for 0 → 4095
// reg pipeline_en;

// always @(posedge clk or negedge rst_n) begin
//     if (!rst_n) begin
//         count       <= 0;
//         pipeline_en <= 0;
//     end else begin

//         if (valid_in) begin
//             // case 1: valid input arrives
//             count       <= 0;
//             pipeline_en <= 1;

//         end else if (count < 4095) begin
//             // case 2: flushing pipeline
//             count       <= count + 1;
//             pipeline_en <= 1;

//         end else if (!valid_in) begin
//             // case 3: finished flushing
//             count       <= 0;
//             pipeline_en <= 0;

//         end

//     end
// end

// module controlunit_4096 #(
//     parameter N = 4096 
// ) (
//     input wire clk,
//     input wire rst_n,
//     input wire valid_in,     // High when DDS starts sending
//     output wire [11:0] addr, 
//     output wire [11:0] sel,
//     output wire pipeline_en, 
//     output wire valid_out    
// );
//     reg [11:0] count,count_valid;

//     // Shift register to delay the valid signal by N-1 cycles
//     reg [N-2:0] valid_pipe; 

//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             count <= 12'd0;
// //           pipeline_en <= 0;
//             valid_pipe <= 0;
//         end else begin
// //            count <= count + 1'b1;
//             if (valid_in) begin
//             // case 1: valid input arrives
//             count_valid       <= 0;
//             count <= count + 1'b1;
// //            pipeline_en <= 1;

//         end else if (count < 4095) begin
//             // case 2: flushing pipeline
//             count_valid       <= count_valid + 1;
// //            pipeline_en <= 1;

//         end else if (!valid_in) begin
//             // case 3: finished flushing
//             count_valid       <= 0;
//             count <= 0;
// //            pipeline_en <= 0;

//         end
//             // Standard SDF delay: Output valid is input valid delayed by N-1
//             valid_pipe <= {valid_pipe[N-3:0], valid_in}; 
//         end
//     end

//     assign sel = count;
//     assign addr = count;
// //    assign pipeline_en = (valid_in || ((count<4095)&&(count>0)))?1:0; // Always running to allow overlap
//     assign pipeline_en = 1;
//     assign valid_out = valid_pipe[N-2]; 
// endmodule/

// module controlunit_4096 #(
//     parameter N = 4096 
// ) (
//     input wire clk,
//     input wire rst_n,
//     input wire valid_in,     // High when DDS starts sending
//     output wire [11:0] addr, 
//     output wire [11:0] sel,
//     output wire pipeline_en, 
//     output wire valid_out
// );
//     reg [11:0] count;
//     // Shift register to delay the valid signal by N-1 cycles
//     reg [N-2:0] valid_pipe; 

//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             count <= 12'd0;
//             valid_pipe <= 0;
//         end else begin
//             count <= count + 1'b1; 
//             // Standard SDF delay: Output valid is input valid delayed by N-1
//             valid_pipe <= {valid_pipe[N-3:0], valid_in}; 
//         end
//     end

//     assign sel = count;
//     assign addr = count;
//     assign pipeline_en = 1'b1; // Always running to allow overlap
//     assign valid_out = valid_pipe[N-2]; 
// endmodule

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