// Analog Devices 
// GP Ain-shams University
// Butterfly unit for 2048-FFT 


/* Description ..........
    This module now drives the pipeline by counting every clock cycle.
*/

// module controlunit_2048 #(
//     parameter N = 2048 
// ) (
//     input  wire        clk,
//     input  wire        rst_n,
//     input  wire        valid_in,
//     output wire [10:0] addr, 
//     output wire [10:0] sel,
//     output wire        pipeline_en, 
//     output wire        valid_out
// );

//     // --------------------------------------------------
//     // Internal Signals
//     // --------------------------------------------------
//     reg [10:0] count;
//     reg [N-2:0] valid_pipe;

//     // --------------------------------------------------
//     // Output valid (delay = N-1 cycles)
//     // --------------------------------------------------
//     assign valid_out = valid_pipe[N-2];

//     // --------------------------------------------------
//     // Pipeline Enable (active when data flows)
//     // --------------------------------------------------
//     assign pipeline_en = valid_in | valid_out;

//     // --------------------------------------------------
//     // Sequential Logic
//     // --------------------------------------------------
//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             count      <= 11'd0;
//             valid_pipe <= { (N-1){1'b0} };
//         end else begin
//             // Shift valid signal through pipeline
//             valid_pipe <= {valid_pipe[N-3:0], valid_in};

//             // Counter follows data movement
//             if (pipeline_en)
//                 count <= count + 1'b1;
//             else
//                 count <= 11'd0;   // reset when idle (frame alignment)
//         end
//     end

//     // --------------------------------------------------
//     // Outputs
//     // --------------------------------------------------
//     assign addr = count +1;
//     assign sel  = count;

// endmodule




// module controlunit_2048 #(
//     parameter N = 2048
// )(
//     input  wire        clk,
//     input  wire        rst_n,
//     input  wire        valid_in,

//     output wire [10:0] sel,
//     output wire [120:0] addr_bus,

//     output wire        pipeline_en,
//     output reg         valid_out
// );

//     // --------------------------------------------------
//     // Internal Signals
//     // --------------------------------------------------
//     reg [10:0] count;
//     reg [10:0] count_pipe [0:9];
//     integer i;

//     // --------------------------------------------------
//     // 1. FLUSH TIMER
//     // --------------------------------------------------
//     // 2048 FFT + 10 pipeline stages - 1
//     // = 2057 cycles flush
//     reg [11:0] flush_counter;

//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             flush_counter <= 0;
//         end
//         else if (valid_in) begin
//             flush_counter <= 12'd2057;
//         end
//         else if (flush_counter > 0) begin
//             flush_counter <= flush_counter - 1'b1;
//         end
//     end

//     // Pipeline active during input OR flushing
//     assign pipeline_en = valid_in | (flush_counter > 0);

//     // --------------------------------------------------
//     // 2. TRAVELING VALID FLAG
//     // --------------------------------------------------
//     reg [2057:0] valid_shift;

//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             valid_shift <= 0;
//             valid_out   <= 0;
//         end
//         else if (pipeline_en) begin
//             valid_shift <= {valid_shift[2056:0], valid_in};
//             valid_out   <= valid_shift[2057];
//         end
//         else begin
//             valid_out <= 0;
//         end
//     end

//     // --------------------------------------------------
//     // 3. COUNTERS
//     // --------------------------------------------------
//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             count <= 0;

//             for (i = 0; i < 10; i = i + 1)
//                 count_pipe[i] <= 0;
//         end
//         else if (pipeline_en) begin
//             count <= count + 1'b1;

//             count_pipe[0] <= count;

//             for (i = 1; i < 10; i = i + 1)
//                 count_pipe[i] <= count_pipe[i-1];
//         end
//     end

//     // --------------------------------------------------
//     // 4. ADVANCED ADDRESS BUS
//     // --------------------------------------------------

//     // Stage 1
//     wire [10:0] p0_plus_1 = count + 1'b1;
//     assign sel[10] = count[10];
//     assign addr_bus[10:0] = p0_plus_1;

//     // Stage 2
//     wire [10:0] p1_plus_1 = count_pipe[0] + 1'b1;
//     assign sel[9] = count_pipe[0][9];
//     assign addr_bus[21:11] = p1_plus_1;

//     // Stage 3
//     wire [10:0] p2_plus_1 = count_pipe[1] + 1'b1;
//     assign sel[8] = count_pipe[1][8];
//     assign addr_bus[32:22] = p2_plus_1;

//     // Stage 4
//     wire [10:0] p3_plus_1 = count_pipe[2] + 1'b1;
//     assign sel[7] = count_pipe[2][7];
//     assign addr_bus[43:33] = p3_plus_1;

//     // Stage 5
//     wire [10:0] p4_plus_1 = count_pipe[3] + 1'b1;
//     assign sel[6] = count_pipe[3][6];
//     assign addr_bus[54:44] = p4_plus_1;

//     // Stage 6
//     wire [10:0] p5_plus_1 = count_pipe[4] + 1'b1;
//     assign sel[5] = count_pipe[4][5];
//     assign addr_bus[65:55] = p5_plus_1;

//     // Stage 7
//     wire [10:0] p6_plus_1 = count_pipe[5] + 1'b1;
//     assign sel[4] = count_pipe[5][4];
//     assign addr_bus[76:66] = p6_plus_1;

//     // Stage 8
//     wire [10:0] p7_plus_1 = count_pipe[6] + 1'b1;
//     assign sel[3] = count_pipe[6][3];
//     assign addr_bus[87:77] = p7_plus_1;

//     // Stage 9
//     wire [10:0] p8_plus_1 = count_pipe[7] + 1'b1;
//     assign sel[2] = count_pipe[7][2];
//     assign addr_bus[98:88] = p8_plus_1;

//     // Stage 10
//     wire [10:0] p9_plus_1 = count_pipe[8] + 1'b1;
//     assign sel[1] = count_pipe[8][1];
//     assign addr_bus[109:99] = p9_plus_1;

//     // Stage 11
//     wire [10:0] p10_plus_1 = count_pipe[9] + 1'b1;
//     assign sel[0] = count_pipe[9][0];
//     assign addr_bus[120:110] = p10_plus_1;

// endmodule


module controlunit_2048 #(
    parameter N = 2048
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in,
    output wire [10:0] sel,             
    output wire [120:0] addr_bus,       
    output wire pipeline_en,
    output reg valid_out
);
    reg [10:0] count;
    (* dont_touch = "true" *) reg [10:0] count_pipe [0:9]; 
    integer i;

    // --- 1. THE FLUSH TIMER ---
    // We need to keep the pipeline running for exactly 2058 cycles
    // AFTER valid_in goes low to push the last piece of data out.
    reg [11:0] flush_counter;
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            flush_counter <= 0;
        end else if (valid_in) begin
            flush_counter <= 12'd2058; // Refresh timer as long as data arrives
        end else if (flush_counter > 0) begin
            flush_counter <= flush_counter - 1; // Count down when valid_in drops
        end
    end

    // Pipeline runs if data is coming IN, or if we are FLUSHING it OUT.
    assign pipeline_en = valid_in | (flush_counter > 0);

    // --- 2. TRAVELING VALID FLAG ---
    // A shift register perfectly delays the valid signal by 2058 cycles.
    // It only shifts when the pipeline is actively running.
    reg [2057:0] valid_shift;
    always @(posedge clk) begin
        // if (!rst_n) begin
        //     valid_shift <= 0;
        //     valid_out <= 0;
        // end else 
        if (pipeline_en) begin
            valid_shift <= {valid_shift[2056:0], valid_in};
            valid_out <= valid_shift[2057]; // Output the delayed valid signal -------------------------------------
        end
        else           valid_out <= 0; // Ensure valid_out is low when pipeline is idle
    end

    // --- 3. COUNTERS ---
    always @(posedge clk) begin
        if (!rst_n) begin
            count <= 0;
            for (i = 0; i < 10; i = i + 1)
                count_pipe[i] <= 0;
        end else if (pipeline_en) begin
            count <= count + 1;
            count_pipe[0] <= count;
            for (i = 1; i < 10; i = i + 1)
                count_pipe[i] <= count_pipe[i-1];
        end
    end

    // --- 4. ADVANCED ADDRESS BUS ---
    // The address sent to the ROMs is advanced by +1 cycle
    // This compensates for the 1-cycle register delay inside the Twiddle ROMs
    wire [10:0] p0_plus_1 = count + 1;
    assign sel[10] = count[10];
    assign addr_bus[10:0] = p0_plus_1;

    wire [10:0] p1_plus_1 = count_pipe[0] + 1;
    assign sel[9] = count_pipe[0][9];
    assign addr_bus[21:11] = p1_plus_1;

    wire [10:0] p2_plus_1 = count_pipe[1] + 1;
    assign sel[8] = count_pipe[1][8];
    assign addr_bus[32:22] = p2_plus_1;

    wire [10:0] p3_plus_1 = count_pipe[2] + 1;
    assign sel[7] = count_pipe[2][7];
    assign addr_bus[43:33] = p3_plus_1;

    wire [10:0] p4_plus_1 = count_pipe[3] + 1;
    assign sel[6] = count_pipe[3][6];
    assign addr_bus[54:44] = p4_plus_1;

    wire [10:0] p5_plus_1 = count_pipe[4] + 1;
    assign sel[5] = count_pipe[4][5];
    assign addr_bus[65:55] = p5_plus_1;

    wire [10:0] p6_plus_1 = count_pipe[5] + 1;
    assign sel[4] = count_pipe[5][4];
    assign addr_bus[76:66] = p6_plus_1;

    wire [10:0] p7_plus_1 = count_pipe[6] + 1;
    assign sel[3] = count_pipe[6][3];
    assign addr_bus[87:77] = p7_plus_1;

    wire [10:0] p8_plus_1 = count_pipe[7] + 1;
    assign sel[2] = count_pipe[7][2];
    assign addr_bus[98:88] = p8_plus_1;

    wire [10:0] p9_plus_1 = count_pipe[8] + 1;
    assign sel[1] = count_pipe[8][1];
    assign addr_bus[109:99] = p9_plus_1;

    wire [10:0] p10_plus_1 = count_pipe[9] + 1;
    assign sel[0] = count_pipe[9][0];
    assign addr_bus[120:110] = p10_plus_1;

endmodule