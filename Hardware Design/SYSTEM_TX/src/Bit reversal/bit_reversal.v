// module bit_reversal_pingpong #(
//     parameter WL = 16,
//     parameter N = 4096,
//     parameter ADDR_W = 12
// )(
//     input clk,
//     input rst_n,
//     input valid_in,
//     input signed [WL-1:0] in_real,
//     input signed [WL-1:0] in_imag,

//     output reg valid_out,
//     output reg signed [WL-1:0] out_real,
//     output reg signed [WL-1:0] out_imag
// );

//     // =========================
//     // Memory (Ping-Pong Banks)
//     // =========================
//     reg signed [WL-1:0] ram_ping_re [0:N-1];
//     reg signed [WL-1:0] ram_ping_im [0:N-1];
//     reg signed [WL-1:0] ram_pong_re [0:N-1];
//     reg signed [WL-1:0] ram_pong_im [0:N-1];

//     // =========================
//     // Control Signals
//     // =========================
//     reg [ADDR_W-1:0] wr_ptr, rd_ptr;
//     reg bank_sel;     // write bank selector
//     reg rd_bank;      // latched read bank
//     reg reading;      // read enable (frame active)

//     // =========================
//     // Bit Reversal Function
//     // =========================
//     function [ADDR_W-1:0] rev;
//         input [ADDR_W-1:0] a;
//         integer i;
//         begin
//             for (i = 0; i < ADDR_W; i = i + 1)
//                 rev[ADDR_W-1-i] = a[i];
//         end
//     endfunction

//     // =========================
//     // Control Logic
//     // =========================
//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             wr_ptr   <= 0;
//             rd_ptr   <= 0;
//             bank_sel <= 0;
//             rd_bank  <= 0;
//             reading  <= 0;
//         end else begin
//             // ================= WRITE =================
//             if (valid_in) begin
//                 wr_ptr <= wr_ptr + 1;

//                 if (wr_ptr == N-1) begin
//                     wr_ptr   <= 0;

//                     // latch which bank to read from
//                     rd_bank  <= bank_sel;

//                     // toggle write bank
//                     bank_sel <= ~bank_sel;

//                     // start reading new frame
//                     reading  <= 1;
//                     rd_ptr   <= 0;
//                 end
//             end

//             // ================= READ =================
//             if (reading) begin
//                 rd_ptr <= rd_ptr + 1;

//                 if (rd_ptr == N-1) begin
//                     reading <= 0;  // stop after full frame
//                 end
//             end
//         end
//     end

//     // =========================
//     // Memory Write
//     // =========================
//     always @(posedge clk) begin
//         if (valid_in) begin
//             if (bank_sel == 0) begin
//                 ram_ping_re[wr_ptr] <= in_real;
//                 ram_ping_im[wr_ptr] <= in_imag;
//             end else begin
//                 ram_pong_re[wr_ptr] <= in_real;
//                 ram_pong_im[wr_ptr] <= in_imag;
//             end
//         end
//     end

//     // =========================
//     // Memory Read
//     // =========================
//     wire [ADDR_W-1:0] rd_addr;
//     assign rd_addr = rev(rd_ptr);

//     always @(posedge clk) begin
//         if (rd_bank == 0) begin
//             out_real <= ram_ping_re[rd_addr];
//             out_imag <= ram_ping_im[rd_addr];
//         end else begin
//             out_real <= ram_pong_re[rd_addr];
//             out_imag <= ram_pong_im[rd_addr];
//         end
//     end

//     // =========================
//     // Valid Output (Aligned)
//     // =========================
//     reg reading_d;

//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             reading_d <= 0;
//             valid_out <= 0;
//         end else begin
//             reading_d <= reading;   // pipeline delay (RAM latency)
//             valid_out <= reading_d;
//         end
//     end

// endmodule
module bit_reversal_pingpong #(
    parameter WL = 16,
    parameter N = 4096,
    parameter ADDR_W = 12
)(
    input clk,
    input rst_n,
    input valid_in,
    input signed [WL-1:0] in_real,
    input signed [WL-1:0] in_imag,

    output reg valid_out,
    output reg signed [WL-1:0] out_real,
    output reg signed [WL-1:0] out_imag
);

    // =========================
    // Memory (Ping-Pong)
    // =========================
    reg signed [WL-1:0] ram_ping_re [0:N-1];
    reg signed [WL-1:0] ram_ping_im [0:N-1];
    reg signed [WL-1:0] ram_pong_re [0:N-1];
    reg signed [WL-1:0] ram_pong_im [0:N-1];

    // =========================
    // Control
    // =========================
    reg [ADDR_W-1:0] wr_ptr, rd_ptr;
    reg bank_sel, rd_bank;
    reg reading;

    // =========================
    // Bit Reversal Function
    // =========================
    function [ADDR_W-1:0] rev;
        input [ADDR_W-1:0] a;
        integer i;
        begin
            for (i = 0; i < ADDR_W; i = i + 1)
                rev[ADDR_W-1-i] = a[i];
        end
    endfunction

    // =========================
    // Control Logic
    // =========================
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            wr_ptr   <= 0;
            rd_ptr   <= 0;
            bank_sel <= 0;
            rd_bank  <= 0;
            reading  <= 0;
        end else begin
            // WRITE
            if (valid_in) begin
                wr_ptr <= wr_ptr + 1;

                if (wr_ptr == N-1) begin
                    wr_ptr   <= 0;
                    rd_bank  <= bank_sel;   // latch read bank
                    bank_sel <= ~bank_sel;  // switch write bank

                    reading  <= 1;
                    rd_ptr   <= 0;
                end
            end

            // READ
            if (reading) begin
                rd_ptr <= rd_ptr + 1;

                if (rd_ptr == N-1)
                    reading <= 0;
            end
        end
    end

    // =========================
    // WRITE
    // =========================
    always @(posedge clk) begin
        if (valid_in) begin
            if (bank_sel == 0) begin
                ram_ping_re[wr_ptr] <= in_real;
                ram_ping_im[wr_ptr] <= in_imag;
            end else begin
                ram_pong_re[wr_ptr] <= in_real;
                ram_pong_im[wr_ptr] <= in_imag;
            end
        end
    end

    // =========================
    // READ (1-cycle latency)
    // =========================
    wire [ADDR_W-1:0] rd_addr = rev(rd_ptr);

    always @(posedge clk) begin
        if (rd_bank == 0) begin
            out_real <= ram_ping_re[rd_addr];
            out_imag <= ram_ping_im[rd_addr];
        end else begin
            out_real <= ram_pong_re[rd_addr];
            out_imag <= ram_pong_im[rd_addr];
        end
    end

    // =========================
    // VALID ALIGNMENT (1 cycle)
    // =========================
 //   reg reading_d;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
//            reading_d <= 0;
            valid_out <= 0;
        end else begin
//            reading_d <= reading;
   //         valid_out <= reading_d;
   valid_out <= reading;
        end
    end

endmodule