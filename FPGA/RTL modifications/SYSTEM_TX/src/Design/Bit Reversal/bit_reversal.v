module bit_reversal_pingpong #(
    parameter WL = 16,
    parameter N = 4096,
    parameter ADDR_W = 12
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in,
    input wire signed [WL-1:0] in_real,
    input wire signed [WL-1:0] in_imag,

    output reg valid_out,
    output wire signed [WL-1:0] out_real,  // Changed to wire to map from BRAM
    output wire signed [WL-1:0] out_imag   // Changed to wire to map from BRAM
);

    // =========================
    // Control Registers
    // =========================
    reg [ADDR_W-1:0] wr_ptr, rd_ptr;
    reg bank_sel, rd_bank, rd_bank_d1;
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
            wr_ptr     <= 0;
            rd_ptr     <= 0;
            bank_sel   <= 0;
            rd_bank    <= 0;
            rd_bank_d1 <= 0;
            reading    <= 0;
            valid_out  <= 0;
        end else begin
            // Shift rd_bank to align with 1-cycle BRAM latency
            rd_bank_d1 <= rd_bank;
            
            // Align valid_out with the 1-cycle BRAM latency
            valid_out <= reading;

            // ---------------------
            // WRITE LOGIC
            // ---------------------
            if (valid_in) begin
                if (wr_ptr == N-1) begin
                    wr_ptr   <= 0;
                    rd_bank  <= bank_sel;   // Latch the bank we just finished writing to
                    bank_sel <= ~bank_sel;  // Switch the write bank
                    reading  <= 1;          // Trigger the read process
                end else begin
                    wr_ptr <= wr_ptr + 1;
                end
            end

            // ---------------------
            // READ LOGIC
            // ---------------------
            if (reading) begin
                if (rd_ptr == N-1) begin
                    rd_ptr <= 0;
                    // Only stop reading if a new frame ISN'T instantly starting
                    if (!(valid_in && wr_ptr == N-1)) begin
                        reading <= 0;
                    end
                end else begin
                    rd_ptr <= rd_ptr + 1;
                end
            end
        end
    end

    // =========================
    // Vivado BRAM Interfacing
    // =========================
    wire [ADDR_W-1:0] rd_addr = rev(rd_ptr);
    
    // Concatenate Real and Imaginary (16 + 16 = 32 bits)
    wire [(WL*2)-1:0] wr_data = {in_real, in_imag};
    wire [(WL*2)-1:0] ping_rd_data, pong_rd_data;
    
    // Write Enables (Port A)
    wire ping_we = valid_in & (bank_sel == 1'b0);
    wire pong_we = valid_in & (bank_sel == 1'b1);

    // -------------------------
    // Instantiate Ping BRAM IP
    // -------------------------
    bram_sdp_32x4096 ping_bram (
        .clka  (clk),
        .wea   (ping_we),
        .addra (wr_ptr),
        .dina  (wr_data),
        
        .clkb  (clk),
        .addrb (rd_addr),
        .doutb (ping_rd_data)
    );

    // -------------------------
    // Instantiate Pong BRAM IP
    // -------------------------
    bram_sdp_32x4096 pong_bram (
        .clka  (clk),
        .wea   (pong_we),
        .addra (wr_ptr),
        .dina  (wr_data),
        
        .clkb  (clk),
        .addrb (rd_addr),
        .doutb (pong_rd_data)
    );

    // =========================
    // Output Multiplexer
    // =========================
    // Uses rd_bank_d1 to wait for the BRAM data to physically exit the IP
    wire [(WL*2)-1:0] read_data_mux = (rd_bank_d1 == 1'b0) ? ping_rd_data : pong_rd_data;

    // Split the 32-bit word back into Real and Imaginary outputs
    assign out_real = read_data_mux[(WL*2)-1 : WL];
    assign out_imag = read_data_mux[WL-1 : 0];

endmodule