module bit_reversal_pingpong #(parameter WL=16, N=4096, ADDR_W=12) (
    input clk, rst_n, valid_in,
    input signed [WL-1:0] in_real, in_imag,
    output reg valid_out,
    output reg signed [WL-1:0] out_real, out_imag
);
    // Dual Memories
    reg signed [WL-1:0] ram_ping_re [0:N-1], ram_ping_im [0:N-1];
    reg signed [WL-1:0] ram_pong_re [0:N-1], ram_pong_im [0:N-1];

    reg [ADDR_W-1:0] wr_ptr, rd_ptr;
    reg bank_sel; // 0: Write Ping/Read Pong, 1: Write Pong/Read Ping
    reg active_read;

    // Bit reversal function stays the same
    function [ADDR_W-1:0] rev;
        input [ADDR_W-1:0] a;
        integer i;
        begin
            for(i=0; i<ADDR_W; i=i+1) rev[ADDR_W-1-i] = a[i];
        end
    endfunction

    always @(posedge clk or negedge rst_n) begin
        if(!rst_n) begin
            wr_ptr <= 0; rd_ptr <= 0; bank_sel <= 0; active_read <= 0;
        end else begin
            if(valid_in) begin
                wr_ptr <= wr_ptr + 1;
                if(wr_ptr == N-1) begin
                    bank_sel <= ~bank_sel;
                    active_read <= 1; // Start reading once the first bank is full
                end
            end
            if(active_read) rd_ptr <= rd_ptr + 1;
        end
    end

    // Memory Logic
    always @(posedge clk) begin
        if(valid_in) begin
            if(!bank_sel) begin ram_ping_re[wr_ptr] <= in_real; ram_ping_im[wr_ptr] <= in_imag; end
            else          begin ram_pong_re[wr_ptr] <= in_real; ram_pong_im[wr_ptr] <= in_imag; end
        end
        valid_out <= active_read;
        if(!bank_sel) begin out_real <= ram_pong_re[rev(rd_ptr)]; out_imag <= ram_pong_im[rev(rd_ptr)]; end
        else          begin out_real <= ram_ping_re[rev(rd_ptr)]; out_imag <= ram_ping_im[rev(rd_ptr)]; end
    end
endmodule