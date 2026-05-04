module controlunit_4096 #(
    parameter N = 4096
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in,
    output wire [11:0] sel,             
    output wire [143:0] addr_bus,       
    output wire pipeline_en,
    output reg valid_out
);
    reg [11:0] count;
    reg [11:0] count_pipe [0:10]; 
    integer i;

    // --- 1. THE FLUSH TIMER ---
    // We need to keep the pipeline running for exactly 4107 cycles
    // AFTER valid_in goes low to push the last piece of data out.
    reg [12:0] flush_counter;
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            flush_counter <= 0;
        end else if (valid_in) begin
            flush_counter <= 13'd4107; // Refresh timer as long as data arrives
        end else if (flush_counter > 0) begin
            flush_counter <= flush_counter - 1; // Count down when valid_in drops
        end
    end

    // Pipeline runs if data is coming IN, or if we are FLUSHING it OUT.
    assign pipeline_en = valid_in | (flush_counter > 0);

    // --- 2. TRAVELING VALID FLAG ---
    // A shift register perfectly delays the valid signal by 4107 cycles.
    // It only shifts when the pipeline is actively running.
    reg [4106:0] valid_shift;
    always @(posedge clk) begin
        // if (!rst_n) begin
        //     valid_shift <= 0;
        // end else 
        if (pipeline_en) begin
            valid_shift <= {valid_shift[4105:0], valid_in};
            valid_out = valid_shift[4106];
        end
    end

    // --- 3. COUNTERS ---
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            count <= 0;
            for (i = 0; i < 11; i = i + 1)
                count_pipe[i] <= 0;
        end else if (pipeline_en) begin
            count <= count + 1;
            count_pipe[0] <= count;
            for (i = 1; i < 11; i = i + 1)
                count_pipe[i] <= count_pipe[i-1];
        end
    end

    // --- 4. ADVANCED ADDRESS BUS ---
    // The address sent to the ROMs is advanced by +1 cycle
    // This compensates for the 1-cycle register delay inside the Twiddle ROMs
    wire [11:0] p0_plus_1 = count + 1;
    assign sel[11] = count[11];
    assign addr_bus[11:0] = p0_plus_1;

    wire [11:0] p1_plus_1 = count_pipe[0] + 1;
    assign sel[10] = count_pipe[0][10];
    assign addr_bus[23:12] = p1_plus_1;

    wire [11:0] p2_plus_1 = count_pipe[1] + 1;
    assign sel[9] = count_pipe[1][9];
    assign addr_bus[35:24] = p2_plus_1;

    wire [11:0] p3_plus_1 = count_pipe[2] + 1;
    assign sel[8] = count_pipe[2][8];
    assign addr_bus[47:36] = p3_plus_1;

    wire [11:0] p4_plus_1 = count_pipe[3] + 1;
    assign sel[7] = count_pipe[3][7];
    assign addr_bus[59:48] = p4_plus_1;

    wire [11:0] p5_plus_1 = count_pipe[4] + 1;
    assign sel[6] = count_pipe[4][6];
    assign addr_bus[71:60] = p5_plus_1;

    wire [11:0] p6_plus_1 = count_pipe[5] + 1;
    assign sel[5] = count_pipe[5][5];
    assign addr_bus[83:72] = p6_plus_1;

    wire [11:0] p7_plus_1 = count_pipe[6] + 1;
    assign sel[4] = count_pipe[6][4];
    assign addr_bus[95:84] = p7_plus_1;

    wire [11:0] p8_plus_1 = count_pipe[7] + 1;
    assign sel[3] = count_pipe[7][3];
    assign addr_bus[107:96] = p8_plus_1;

    wire [11:0] p9_plus_1 = count_pipe[8] + 1;
    assign sel[2] = count_pipe[8][2];
    assign addr_bus[119:108] = p9_plus_1;

    wire [11:0] p10_plus_1 = count_pipe[9] + 1;
    assign sel[1] = count_pipe[9][1];
    assign addr_bus[131:120] = p10_plus_1;

    wire [11:0] p11_plus_1 = count_pipe[10] + 1;
    assign sel[0] = count_pipe[10][0];
    assign addr_bus[143:132] = p11_plus_1;

endmodule