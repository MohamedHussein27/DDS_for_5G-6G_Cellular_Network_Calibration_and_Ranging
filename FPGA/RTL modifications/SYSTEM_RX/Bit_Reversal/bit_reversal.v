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
//     // Memory (Ping-Pong)
//     // =========================
//     reg signed [WL-1:0] ram_ping_re [0:N-1];
//     reg signed [WL-1:0] ram_ping_im [0:N-1];
//     reg signed [WL-1:0] ram_pong_re [0:N-1];
//     reg signed [WL-1:0] ram_pong_im [0:N-1];

//     // =========================
//     // Control
//     // =========================
//     reg [ADDR_W-1:0] wr_ptr, rd_ptr;
//     reg bank_sel, rd_bank;
//     reg reading;

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
// //            valid_out <= 0;
// //            out_real <= 0;
// //            out_imag <= 0;
// //            ram_ping_re[0] <= 0;
// //            ram_ping_im[0] <= 0;
//         end else begin
//             // WRITE
//             if (valid_in) begin
//                 wr_ptr <= wr_ptr + 1;

//                 if (wr_ptr == N-1) begin
//                     wr_ptr   <= 0;
//                     rd_bank  <= bank_sel;   // latch read bank
//                     bank_sel <= ~bank_sel;  // switch write bank

//                     reading  <= 1;
//                     rd_ptr   <= 0;
//                 end
//             end

//             // READ
//             if (reading) begin
//                 rd_ptr <= rd_ptr + 1;

//                 if (rd_ptr == N-1)
//                     reading <= 0;
//             end
//         end
//     end

//     // =========================
//     // WRITE
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
//     // READ
//     // =========================
//     wire [ADDR_W-1:0] rd_addr = rev(rd_ptr);

//     always @(posedge clk) begin
//         if (rd_bank == 0) begin
//             out_real <= ram_ping_re[rd_addr];
//             out_imag <= ram_ping_im[rd_addr];
//         end else begin
//             out_real <= ram_pong_re[rd_addr];
//             out_imag <= ram_pong_im[rd_addr];
//         end
//     end



//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             valid_out <= 0;
//         end else begin
//             valid_out <= reading;
//         end
//     end

// endmodule

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
//     // Memory (Ping-Pong)
//     // =========================
//     // reg signed [WL-1:0] ram_ping_re [0:N-1];
//     // reg signed [WL-1:0] ram_ping_im [0:N-1];
//     // reg signed [WL-1:0] ram_pong_re [0:N-1];
//     // reg signed [WL-1:0] ram_pong_im [0:N-1];

//     (* ram_style = "block" *) reg [WL*2-1:0] ram_ping [0:N-1];
//     (* ram_style = "block" *) reg [WL*2-1:0] ram_pong [0:N-1];

//     // =========================
//     // Control
//     // =========================
//     reg [ADDR_W-1:0] wr_ptr, rd_ptr;
//     reg bank_sel, rd_bank;
//     reg reading;

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
// //            valid_out <= 0;
// //            out_real <= 0;
// //            out_imag <= 0;
// //            ram_ping_re[0] <= 0;
// //            ram_ping_im[0] <= 0;
//         end else begin
//             // WRITE
//             if (valid_in) begin
//                 wr_ptr <= wr_ptr + 1;

//                 if (wr_ptr == N-1) begin
//                     wr_ptr   <= 0;
//                     rd_bank  <= bank_sel;   // latch read bank
//                     bank_sel <= ~bank_sel;  // switch write bank

//                     reading  <= 1;
//                     rd_ptr   <= 0;
//                 end
//             end

//             // READ
//             if (reading) begin
//                 rd_ptr <= rd_ptr + 1;

//                 if (rd_ptr == N-1)
//                     reading <= 0;
//             end
//         end
//     end

//     // // =========================
//     // // WRITE
//     // // =========================
//     // always @(posedge clk) begin
//     //     if (valid_in) begin
//     //         if (bank_sel == 0) begin
//     //             ram_ping_re[wr_ptr] <= in_real;
//     //             ram_ping_im[wr_ptr] <= in_imag;
//     //         end else begin
//     //             ram_pong_re[wr_ptr] <= in_real;
//     //             ram_pong_im[wr_ptr] <= in_imag;
//     //         end
//     //     end
//     // end

//     // // =========================
//     // // READ
//     // // =========================
//     // wire [ADDR_W-1:0] rd_addr = rev(rd_ptr);

//     // always @(posedge clk) begin
//     //     if (rd_bank == 0) begin
//     //         out_real <= ram_ping_re[rd_addr];
//     //         out_imag <= ram_ping_im[rd_addr];
//     //     end else begin
//     //         out_real <= ram_pong_re[rd_addr];
//     //         out_imag <= ram_pong_im[rd_addr];
//     //     end
//     // end

//     // =========================
//     // WRITE
//     // =========================
//     always @(posedge clk) begin
//         if (valid_in) begin
//             if (bank_sel == 0) begin
//                 ram_ping[wr_ptr] <= {in_real, in_imag}; // Combine into 32-bit word
//             end else begin
//                 ram_pong[wr_ptr] <= {in_real, in_imag}; // Combine into 32-bit word
//             end
//         end
//     end

//     // =========================
//     // READ
//     // =========================
//     wire [ADDR_W-1:0] rd_addr = rev(rd_ptr);
//     reg [WL*2-1:0] read_data; // Intermediate 32-bit register

//     always @(posedge clk) begin
//         if (rd_bank == 0) begin
//             read_data <= ram_ping[rd_addr];
//         end else begin
//             read_data <= ram_pong[rd_addr];
//         end
//         // Split the 32-bit word back into Real and Imaginary outputs
//         out_imag <= read_data[WL-1:0];
//         out_real <= read_data[WL*2-1:WL];
//     end



//     always @(posedge clk or negedge rst_n) begin
//         if (!rst_n) begin
//             valid_out <= 0;
//         end else begin
//             valid_out <= reading;
//         end
//     end

// endmodule

`timescale 1ns / 1ps

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
    // COMBINED into 32-bit words (16-bit Imaginary, 16-bit Real)
    // Forced to use BRAM to fix the 17,872 LUTRAM Place Design error
    (* ram_style = "block" *) reg [WL*2-1:0] ram_ping [0:N-1];
    (* ram_style = "block" *) reg [WL*2-1:0] ram_pong [0:N-1];

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
    always @(posedge clk) begin
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
                // if (rd_ptr == N)
                    reading <= 0;
            end
        end
    end

    // =========================
    // WRITE (Combined 32-bit Word)
    // =========================
    // Create explicit Write Enable signals for perfect BRAM inference
    (* dont_touch = "true" *) wire ping_we = valid_in & (bank_sel == 0);
    (* dont_touch = "true" *) wire pong_we = valid_in & (bank_sel == 1);

    always @(posedge clk) begin
        if (ping_we) ram_ping[wr_ptr] <= {in_imag, in_real};
        if (pong_we) ram_pong[wr_ptr] <= {in_imag, in_real};
    end

    // =========================
    // READ (Combined 32-bit Word)
    // =========================
    wire [ADDR_W-1:0] rd_addr = rev(rd_ptr);
    (* dont_touch = "true" *) reg [WL*2-1:0] ping_out;
    (* dont_touch = "true" *) reg [WL*2-1:0] pong_out;

    // Read BOTH memories synchronously into dedicated registers (Infers BRAM DOUT)
    always @(posedge clk) begin
        ping_out <= ram_ping[rd_addr];
        pong_out <= ram_pong[rd_addr];
    end

    // Fabric MUX selects the active bank AFTER the BRAM registers
    wire [WL*2-1:0] read_data = (rd_bank == 0) ? ping_out : pong_out;

    // Symmetrically unpack the word back into 16-bit registers.
    // Using a combinational block here preserves your exact original cycle timing!
    always @(*) begin
        out_real = read_data[WL-1:0];       // Lower 16 bits = Real
        out_imag = read_data[WL*2-1:WL];    // Upper 16 bits = Imaginary
    end

    // =========================
    // VALID OUT
    // =========================
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            valid_out <= 0;
        end else begin
            valid_out <= reading;
        end
    end

endmodule