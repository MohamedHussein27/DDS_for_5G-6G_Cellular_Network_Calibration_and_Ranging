module RX_wrapper #(
    parameter WL = 16,
    parameter N  = 4096
)(
    input  wire                 clk,
    input  wire                 rst_n,
    input  wire                 rx_valid_in,
    output wire                 real_match_flag_ofdm,
    output wire                 imag_match_flag_ofdm,
    output wire                 real_match_flag_radar,
    output wire                 imag_match_flag_radar

);
    // Channel Input

    wire signed [WL-1:0] rx_in_re;
    wire signed [WL-1:0] rx_in_im;
    
    // OFDM Output
    wire                 ofdm_valid_out;
    wire signed [WL-1:0] ofdm_out_re;
    wire signed [WL-1:0] ofdm_out_im;
    
    // Radar Output
    wire                radar_valid_out;
    wire signed [WL-1:0] radar_out_re;
    wire signed [WL-1:0] radar_out_im;

    //REF RAM Outputs for Comparison
    wire signed [WL-1:0] ofdm_ref_re;
    wire signed [WL-1:0] ofdm_ref_im;
    wire signed [WL-1:0] radar_ref_re;
    wire signed [WL-1:0] radar_ref_im;
    // Delayed Write Address for Comparison
    wire [10:0] write_addr_ofdm;
    wire [10:0] write_addr_radar;
    reg [10:0] write_addr_d_ofdm;
    reg [10:0] write_addr_d_radar;
    // Comparison Flags (declared as output ports above; no re-declaration needed)

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            write_addr_d_ofdm <= 0;
            write_addr_d_radar <= 0;
        end else begin
            write_addr_d_ofdm <= write_addr_ofdm; // Delay the address by 1 cycle to match RAM output timing
            write_addr_d_radar <= write_addr_radar; // Delay the address by 1 cycle to match RAM output timing
        end
    end

    reg rx_valid_in_delayed;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) rx_valid_in_delayed <= 0;
        else rx_valid_in_delayed <= rx_valid_in;
    end

    wire [11:0] in_count;

    input_counter #(
        .ADDR_WIDTH(12)
    ) u_input_counter (
        .enable(rx_valid_in),
        .clk(clk),
        .rst_n(rst_n),
        .count(in_count)
    );

    rx_input_rom_re #(
        .WL(WL),
        .DEPTH(N)
    ) u_rx_input_gen_re (
        .clk(clk),
        .addr(in_count),
        .data_out_re(rx_in_re)
    );

    rx_input_rom_im #(
        .WL(WL),
        .DEPTH(N)
    ) u_rx_input_gen_im (
        .clk(clk),
        .addr(in_count),
        .data_out_im(rx_in_im)
    );

    // RX_input_gen_re u_rx_input_gen_re (
    //     .clka(clk),
    //     .addra(in_count),
    //     .douta(rx_in_re)
    // );

    // RX_input_gen_im u_rx_input_gen_im (
    //     .clka(clk),
    //     .addra(in_count),
    //     .douta(rx_in_im)
    // );

    // Instantiate the RX_TOP module
    RX_TOP #(
        .WL(WL),
        .N(N)
    ) u_rx_top (
        .clk(clk),
        .rst_n(rst_n),
        
        // Channel Input
        .rx_valid_in(rx_valid_in_delayed),
        .rx_in_re(rx_in_re),
        .rx_in_im(rx_in_im),
        
        // OFDM Output
        .ofdm_valid_out(ofdm_valid_out),
        .ofdm_out_re(ofdm_out_re),
        .ofdm_out_im(ofdm_out_im),
        
        // Radar Output
        .radar_valid_out(radar_valid_out),
        .radar_out_re(radar_out_re),
        .radar_out_im(radar_out_im)
    );

    bram_write_counter #(
        .ADDR_WIDTH(11)
    ) u_bram_write_counter (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in_ofdm(ofdm_valid_out),
        .valid_in_radar(radar_valid_out),
        .write_addr_ofdm(write_addr_ofdm),
        .write_addr_radar(write_addr_radar)
    );

    OFDM_ref_mem_re u_ofdm_ref_mem_re (
        .clka(clk),
        .addra(write_addr_ofdm),
        .douta(ofdm_ref_re)
    ); 
    OFDM_ref_mem_im u_ofdm_ref_mem_im (
        .clka(clk),
        .addra(write_addr_ofdm),
        .douta(ofdm_ref_im)
    );
    RADAR_ref_mem_re u_radar_ref_mem_re (
        .clka(clk),
        .addra(write_addr_radar),
        .douta(radar_ref_re)
    ); 
    RADAR_ref_mem_im u_radar_ref_mem_im (
        .clka(clk),
        .addra(write_addr_radar),
        .douta(radar_ref_im)
    ); 


    // Feed the DELAYED signals to perfectly match the ROM outputs
    real_comparator #(
        .INPUT_WIDTH(WL), 
        .NUM_SYMBOLS(2048), 
        .SYMBOL(11) 
    ) u_real_comp (
        .clk(clk) , 
        .rst_n(rst_n),
        .ofdm_valid_in(ofdm_valid_out),
        .radar_valid_in(radar_valid_out),
        .OFDM_real_in(ofdm_out_re),
        .OFDM_real_ref_in(ofdm_ref_re),
        .radar_real_in(radar_out_re),
        .radar_real_ref_in(radar_ref_re),
        .symbol_index_ofdm(write_addr_d_ofdm),
        .symbol_index_radar(write_addr_d_radar),
        .comparison_result_ofdm(real_match_flag_ofdm),
        .comparison_result_radar(real_match_flag_radar)
    );

    imaginary_comparator #(
        .INPUT_WIDTH(WL), 
        .NUM_SYMBOLS(2048), 
        .SYMBOL(11)
    ) u_imag_comp (
        .clk(clk),
        .rst_n(rst_n),
        .ofdm_valid_in(ofdm_valid_out),
        .radar_valid_in(radar_valid_out),
        .OFDM_imag_in(ofdm_out_im),
        .OFDM_imag_ref_in(ofdm_ref_im),
        .radar_imag_in(radar_out_im),
        .radar_imag_ref_in(radar_ref_im),
        .symbol_index_ofdm(write_addr_d_ofdm),
        .symbol_index_radar(write_addr_d_radar),
        .comparison_result_ofdm(imag_match_flag_ofdm),
        .comparison_result_radar(imag_match_flag_radar)
    );

endmodule