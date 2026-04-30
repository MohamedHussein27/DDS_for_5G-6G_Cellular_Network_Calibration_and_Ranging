`timescale 1ns / 1ps

module SYSTEM_TOP_tb();

    parameter WL = 16;
    parameter N  = 4096;
    parameter DDS_W = 8;

    reg clk;
    reg rst_n;
    
    reg dds_enable;
    reg [31:0] FTW_start;
    reg [12:0] cycles;
    reg [31:0] FTW_step;

    reg signed [WL-1:0] ofdm_in_re;
    reg signed [WL-1:0] ofdm_in_im;
    wire ofdm_rd_en;

    wire ofdm_valid_out;
    wire signed [WL-1:0] ofdm_out_re;
    wire signed [WL-1:0] ofdm_out_im;
    
    wire radar_valid_out;
    wire signed [WL-1:0] radar_out_re;
    wire signed [WL-1:0] radar_out_im;

    SYSTEM_TOP #(
        .WL(WL), .N(N)
 //       , .DDS_W(DDS_W)
    ) uut (
        .clk(clk), .rst_n(rst_n),
        .dds_enable(dds_enable),
        .FTW_start(FTW_start), .cycles(cycles), .FTW_step(FTW_step),
        .ofdm_in_re(ofdm_in_re), .ofdm_in_im(ofdm_in_im),
        .ofdm_rd_en(ofdm_rd_en),
        .ofdm_valid_out(ofdm_valid_out),
        .ofdm_out_re(ofdm_out_re), .ofdm_out_im(ofdm_out_im),
        .radar_valid_out(radar_valid_out),
        .radar_out_re(radar_out_re), .radar_out_im(radar_out_im)
    );

    // Clock Generation (500 MHz)
    initial begin
        clk = 0;
        forever #1 clk = ~clk; 
    end

    // OFDM Memory Loading & Streaming
    reg signed [WL-1:0] ofdm_rom_re [0:2047];
    reg signed [WL-1:0] ofdm_rom_im [0:2047];
    reg [11:0] ofdm_ptr;

    initial begin
        $readmemh("ofdm_data_re.hex", ofdm_rom_re);
        $readmemh("ofdm_data_im.hex", ofdm_rom_im);
        ofdm_ptr = 0;
    end

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

    // Output File Pointers & Counters
    integer f_rx_ofdm_re, f_rx_ofdm_im;
    integer f_rx_radar_re, f_rx_radar_im;
    integer f_ofdm_out_re, f_ofdm_out_im;
    
    integer ofdm_cnt  = 0;
    integer radar_cnt = 0;
    
    // --- ADD THIS LINE ---
    integer f_mult_re, f_mult_im;

    initial begin
        f_rx_ofdm_re  = $fopen("rtl_rx_ofdm_re.hex", "w");
        f_rx_ofdm_im  = $fopen("rtl_rx_ofdm_im.hex", "w");
        f_rx_radar_re = $fopen("rtl_rx_radar_re.hex", "w");
        f_rx_radar_im = $fopen("rtl_rx_radar_im.hex", "w");
        
        f_ofdm_out_re = $fopen("final_ofdm_re.hex", "w");
        f_ofdm_out_im = $fopen("final_ofdm_im.hex", "w");

        // --- ADD THESE ---
        f_mult_re = $fopen("rtl_mult_re.hex", "w");
        f_mult_im = $fopen("rtl_mult_im.hex", "w");

        if (f_ofdm_out_re == 0 || f_mult_re == 0) begin
            $display("ERROR: Could not open one or more hex files for writing.");
            $finish;
        end

        rst_n      = 0;
        dds_enable = 0;
        FTW_start  = 32'd0; 
        FTW_step   = 32'd426666;
        cycles     = 13'd4096;

        #20 rst_n = 1;
        #20 dds_enable = 1;

        #1000000;
        $display("ERROR: Simulation timed out!");
        $stop;
    end

    // =========================================================
    // LOGGING LOGIC: FINAL OFDM OUTPUT (FREQUENCY DOMAIN)
    // =========================================================
    // Capture processed symbols at the end of the chain
    always @(posedge clk) begin
        if (ofdm_valid_out) begin
            $fdisplay(f_ofdm_out_re, "%04x", ofdm_out_re & 16'hFFFF);
            $fdisplay(f_ofdm_out_im, "%04x", ofdm_out_im & 16'hFFFF);
            
            $fdisplay(f_rx_ofdm_re, "%04x", ofdm_out_re);
            $fdisplay(f_rx_ofdm_im, "%04x", ofdm_out_im);
            ofdm_cnt = ofdm_cnt + 1;
        end
    end

    // =========================================================
    // LOGGING LOGIC: INTERNAL MULTIPLIER OUTPUT
    // =========================================================
    // Probe into the RX_TOP instance to grab the intermediate conjugate multiplication
// Probing into the RX_TOP instance to grab the intermediate conjugate multiplication
// Probing into the RX_TOP instance
    always @(posedge clk) begin
        if (uut.u_rx_top.mult_valid) begin
            $fdisplay(f_mult_re, "%04x", uut.u_rx_top.mult_re & 16'hFFFF);
            $fdisplay(f_mult_im, "%04x", uut.u_rx_top.mult_im & 16'hFFFF);
        end
    end

    // =========================================================
    // LOGGING LOGIC: RX RADAR OUTPUT
    // =========================================================
    // Capture RX Radar Output (4096 Samples)
    always @(posedge clk) begin
        if (radar_valid_out) begin
            $fdisplay(f_rx_radar_re, "%04x", radar_out_re);
            $fdisplay(f_rx_radar_im, "%04x", radar_out_im);
            radar_cnt = radar_cnt + 1;
            
            // The Radar output is the very last block in the pipeline
            if (radar_cnt == 2048) begin
                $display("SUCCESS: Full System Transmit & Receive Complete!");
                $fclose(f_rx_ofdm_re);  $fclose(f_rx_ofdm_im);
                $fclose(f_rx_radar_re); $fclose(f_rx_radar_im);
                $fclose(f_ofdm_out_re); $fclose(f_ofdm_out_im);
                
                // --- NEW: Close multiplier logs ---
                $fclose(f_mult_re);     $fclose(f_mult_im);
                #100; 
                $stop;
            end
        end
    end

endmodule