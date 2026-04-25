`timescale 1ns / 1ps

module TX_TOP_tb();

    // --- Parameters ---
    parameter WL = 16;
    parameter N  = 4096;
    parameter DDS_W = 8;

    // --- Testbench Signals ---
    reg clk;
    reg rst_n;
    
    // DDS Controls
    reg dds_enable;
    reg [31:0] FTW_start;
    reg [12:0] cycles;
    reg [31:0] FTW_step;

    // OFDM Data Interface
    reg signed [WL-1:0] ofdm_in_re;
    reg signed [WL-1:0] ofdm_in_im;
    wire ofdm_rd_en;

    // Output Interface
    wire tx_valid;
    wire signed [WL-1:0] tx_out_re;
    wire signed [WL-1:0] tx_out_im;

    // --- Instantiating the TX_TOP ---
    TX_TOP #(
        .WL(WL),
        .N(N),
        .DDS_W(DDS_W)
    ) uut (
        .clk(clk),
        .rst_n(rst_n),
        .dds_enable(dds_enable),
        .FTW_start(FTW_start),
        .cycles(cycles),
        .FTW_step(FTW_step),
        .ofdm_in_re(ofdm_in_re),
        .ofdm_in_im(ofdm_in_im),
        .ofdm_rd_en(ofdm_rd_en),
        .tx_valid(tx_valid),
        .tx_out_re(tx_out_re),
        .tx_out_im(tx_out_im)
    );

    // --- 1. Clock Generation ---
    initial begin
        clk = 0;
        forever #1 clk = ~clk; // 500 MHz simulation clock (2ns period)
    end

    // --- 2. Memory Buffers & File Loading ---
    reg signed [WL-1:0] ofdm_rom_re [0:2047];
    reg signed [WL-1:0] ofdm_rom_im [0:2047];
    reg [11:0] ofdm_ptr;

    initial begin
        $readmemh("ofdm_data_re.hex", ofdm_rom_re);
        $readmemh("ofdm_data_im.hex", ofdm_rom_im);
        ofdm_ptr = 0;
    end

    // Stream the OFDM data
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            ofdm_ptr   <= 0;
            ofdm_in_re <= 0;
            ofdm_in_im <= 0;
        end else if (ofdm_rd_en) begin
            ofdm_in_re <= ofdm_rom_re[ofdm_ptr];
            ofdm_in_im <= ofdm_rom_im[ofdm_ptr];
            ofdm_ptr   <= ofdm_ptr + 1;
        end
    end

    // --- 3. Output Capturing & Simulation Control ---
    integer file_re, file_im, file_dds;
    integer sample_count;
    reg dds_valid_delayed; // Used to detect the falling edge of DDS valid

    initial begin
        // Open files
        file_re  = $fopen("rtl_tx_out_re.txt", "w");
        file_im  = $fopen("rtl_tx_out_im.txt", "w");
        file_dds = $fopen("rtl_dds_out.txt", "w"); // <--- NEW: DDS file
        
        sample_count = 0;
        dds_valid_delayed = 0;

        // Initialize Inputs
        rst_n      = 0;
        dds_enable = 0;
        FTW_start  = 32'd0; 
        FTW_step   = 32'd426666;
        cycles     = 13'd4096;

        // Apply Reset
        #20 rst_n = 1;

        // Start the Transmitter
        #20 dds_enable = 1;

        // Timeout protection
        #500000;
        $display("ERROR: Simulation timed out!");
        $stop;
    end

    // --- NEW: Capture DDS Internal Output ---
    always @(posedge clk) begin
        if (rst_n) begin
            // Keep track of the previous clock cycle's valid signal
            dds_valid_delayed <= uut.dds_valid;

            // Write to file while valid is HIGH
            if (uut.dds_valid) begin
                // $signed forces Verilog to treat the 8-bit vector as a negative number if the MSB is 1
                $fdisplay(file_dds, "%d", $signed(uut.dds_amplitude));
            end

            // Close file exactly when valid drops from 1 to 0 (Falling Edge)
            if (dds_valid_delayed == 1'b1 && uut.dds_valid == 1'b0) begin
                $display("SUCCESS: Captured DDS samples. Closing DDS file.");
                $fclose(file_dds);
            end
        end
    end

    // --- TX IFFT Final Output Capture ---
    always @(posedge clk) begin
        if (tx_valid) begin
            $fdisplay(file_re, "%d", tx_out_re);
            $fdisplay(file_im, "%d", tx_out_im);
            
            sample_count = sample_count + 1;
            
            if (sample_count == 4096) begin
                $display("SUCCESS: Captured exactly 4096 TX output samples.");
                $fclose(file_re);
                $fclose(file_im);
                #100; // Delay to let the OS save the file
                $stop;
            end
        end
    end

endmodule