`timescale 1ns / 1ps

module TX_TOP_tb();

    // =========================================================================
    // 1. PARAMETERS
    // =========================================================================
    parameter WL = 16;
    parameter N  = 4096;
    parameter DDS_W = 8;

    // =========================================================================
    // 2. TESTBENCH SIGNALS
    // =========================================================================
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

    // =========================================================================
    // 3. DUT INSTANTIATION
    // =========================================================================
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

    // =========================================================================
    // 4. CLOCK GENERATION
    // =========================================================================
    initial begin
        clk = 0;
        forever #1 clk = ~clk; // 500 MHz simulation clock (2ns period)
    end

    // =========================================================================
    // 5. OFDM MEMORY BUFFERING
    // =========================================================================
    reg signed [WL-1:0] ofdm_rom_re [0:2047];
    reg signed [WL-1:0] ofdm_rom_im [0:2047];
    reg [11:0] ofdm_ptr;

    initial begin
        // Make sure these hex files exist in your simulation directory
        $readmemh("ofdm_data_re.hex", ofdm_rom_re);
        $readmemh("ofdm_data_im.hex", ofdm_rom_im);
        ofdm_ptr = 0;
    end

    // Stream the OFDM data when requested by the MUX
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

    // =========================================================================
    // 6. SIMULATION CONTROL & FILE I/O
    // =========================================================================
    integer file_re, file_im;               // Final TX Outputs
    integer file_dds;                       // DDS Time-Domain Output
    integer file_raw_re, file_raw_im;       // Raw Scrambled FFT Output
    integer file_bitrev_re, file_bitrev_im; // Bit Reversal Hex Output
    integer file_mux_re, file_mux_im;       // <--- NEW: MUX Hex Output
    
    integer sample_count;
    integer raw_fft_count;
    integer bitrev_count;
    integer mux_count;                      // <--- NEW: Counter for MUX
    
    reg dds_valid_delayed; 
    reg raw_fft_valid_delayed;
    reg bitrev_valid_delayed;
    reg mux_valid_delayed;                  // <--- NEW: Edge detector for MUX

    initial begin
        // Open all extraction files
        file_re        = $fopen("rtl_tx_out_re.txt", "w");
        file_im        = $fopen("rtl_tx_out_im.txt", "w");
        file_dds       = $fopen("rtl_dds_out.txt", "w"); 
        file_raw_re    = $fopen("rtl_raw_fft_re.txt", "w");
        file_raw_im    = $fopen("rtl_raw_fft_im.txt", "w");
        file_bitrev_re = $fopen("rtl_bitrev_out_re.hex", "w");
        file_bitrev_im = $fopen("rtl_bitrev_out_im.hex", "w");
        
        // --- NEW: Open hex files for MUX extraction ---
        file_mux_re    = $fopen("rtl_mux_out_re.hex", "w");
        file_mux_im    = $fopen("rtl_mux_out_im.hex", "w");
        
        sample_count = 0;
        raw_fft_count = 0;
        bitrev_count = 0;
        mux_count = 0;
        
        dds_valid_delayed = 0;
        raw_fft_valid_delayed = 0;
        bitrev_valid_delayed = 0;
        mux_valid_delayed = 0;

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

    // -------------------------------------------------------------------------
    // CAPTURE 1: Internal DDS Output (Decimal)
    // -------------------------------------------------------------------------
    always @(posedge clk) begin
        if (rst_n) begin
            dds_valid_delayed <= uut.dds_valid;

            if (uut.dds_valid) begin
                $fdisplay(file_dds, "%d", $signed(uut.dds_amplitude));
            end

            if (dds_valid_delayed == 1'b1 && uut.dds_valid == 1'b0) begin
                $display("SUCCESS: Captured DDS samples. Closing DDS file.");
                $fclose(file_dds);
            end
        end
    end

    // -------------------------------------------------------------------------
    // CAPTURE 2: Internal Raw Scrambled FFT Output (Decimal)
    // -------------------------------------------------------------------------
    always @(posedge clk) begin
        if (rst_n) begin
            raw_fft_valid_delayed <= uut.u_fft_tx.valid_out;

            if (uut.u_fft_tx.valid_out) begin
                $fdisplay(file_raw_re, "%d", $signed(uut.u_fft_tx.out_real));
                $fdisplay(file_raw_im, "%d", $signed(uut.u_fft_tx.out_imag));
                raw_fft_count = raw_fft_count + 1;
            end

            if (raw_fft_valid_delayed == 1'b1 && uut.u_fft_tx.valid_out == 1'b0) begin
                $display("SUCCESS: Captured %0d Raw FFT samples. Closing files.", raw_fft_count);
                $fclose(file_raw_re);
                $fclose(file_raw_im);
            end
        end
    end

    // -------------------------------------------------------------------------
    // CAPTURE 3: Internal Bit Reversal Output (HEXADECIMAL)
    // -------------------------------------------------------------------------
    always @(posedge clk) begin
        if (rst_n) begin
            bitrev_valid_delayed <= uut.bit_rev_valid;

            if (uut.bit_rev_valid) begin
                $fdisplay(file_bitrev_re, "%04x", uut.bit_rev_re);
                $fdisplay(file_bitrev_im, "%04x", uut.bit_rev_im);
                bitrev_count = bitrev_count + 1;
            end

            if (bitrev_valid_delayed == 1'b1 && uut.bit_rev_valid == 1'b0) begin
                $display("SUCCESS: Captured %0d Bit Reversal Hex samples. Closing files.", bitrev_count);
                $fclose(file_bitrev_re);
                $fclose(file_bitrev_im);
            end
        end
    end

    // -------------------------------------------------------------------------
    // CAPTURE 4: Internal MUX Output (HEXADECIMAL)
    // -------------------------------------------------------------------------
    always @(posedge clk) begin
        if (rst_n) begin
            // Tap the interconnected wires inside TX_TOP
            mux_valid_delayed <= uut.mux_valid;

            if (uut.mux_valid) begin
                // %04x forces 4 lowercase hex characters, zero-padded (e.g., '11a9')
                $fdisplay(file_mux_re, "%04x", uut.mux_re);
                $fdisplay(file_mux_im, "%04x", uut.mux_im);
                mux_count = mux_count + 1;
            end

            // Close files on the falling edge of mux_valid
            if (mux_valid_delayed == 1'b1 && uut.mux_valid == 1'b0) begin
                $display("SUCCESS: Captured %0d MUX Hex samples. Closing files.", mux_count);
                $fclose(file_mux_re);
                $fclose(file_mux_im);
            end
        end
    end

    // -------------------------------------------------------------------------
    // CAPTURE 5: Final TX Output (Post-IFFT) (Decimal)
    // -------------------------------------------------------------------------
    always @(posedge clk) begin
        if (tx_valid) begin
            $fdisplay(file_re, "%d", tx_out_re);
            $fdisplay(file_im, "%d", tx_out_im);
            
            sample_count = sample_count + 1;
            
            // Stop simulation once the 4096-bin frame is completely transmitted
            if (sample_count == 4096) begin
                $display("SUCCESS: Captured exactly 4096 Final TX output samples.");
                $fclose(file_re);
                $fclose(file_im);
                #100; // Delay to let the OS flush the files
                $stop;
            end
        end
    end

endmodule