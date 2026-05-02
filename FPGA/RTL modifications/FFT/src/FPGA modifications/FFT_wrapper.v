module FFT_wrapper #(
    parameter WL = 16,
    parameter ADDR_WIDTH = 12,
    parameter INPUT_WIDTH = 16,
    parameter NUM_SYMBOLS = 4096,
    parameter SYMBOL = 12
) (
    input wire clk,
    input wire rst_n,
    input wire start,
    //input wire valid_in,
    output wire real_comparison_result,
    output wire imag_comparison_result  
);

    wire [ADDR_WIDTH-1:0] write_addr;
    
    //dds
    wire signed [WL-1:0] dds_douta;

    wire valid_out;
    wire signed [WL-1:0] out_real;
    wire signed [WL-1:0] out_imag;
    wire signed [WL-1:0] fft_imag_ref_in;
    wire signed [WL-1:0] fft_real_ref_in;

    wire [ADDR_WIDTH-1:0] dds_count;

    bram_write_counter #(.ADDR_WIDTH (ADDR_WIDTH)) u_counter (
        .clk (clk),
        .rst_n (rst_n),
        .valid_in (valid_out),
        .write_addr (write_addr)
    );

    DDS_input u_dds (
          .clka (clk),
          .addra (dds_count),
          .douta (dds_douta)
    );

    fft_4096_top #(.WL(WL)) u_fft (
             .clk (clk),
             .rst_n (rst_n),
             .valid_in (start),           //<--- NEW
             .in_real (dds_douta),
             .in_imag (0), 
             .valid_out (valid_out),       //  <--- NEW
             .out_real (out_real),
             .out_imag (out_imag)
    );

    real_comparator #(
     .INPUT_WIDTH (INPUT_WIDTH), .NUM_SYMBOLS (NUM_SYMBOLS), .SYMBOL (SYMBOL) )  u_real_comp (

    .clk (clk) , 
    .rst_n (rst_n),
    .valid_in (valid_out),
    .fft_real_in (out_real),
    .fft_real_ref_in (fft_real_ref_in),
    .symbol_index (write_addr),
    .comparison_result (real_comparison_result)

);
    imaginary_comparator #(.INPUT_WIDTH (INPUT_WIDTH), .NUM_SYMBOLS (NUM_SYMBOLS), .SYMBOL (SYMBOL)) u_imag_comp (

      .clk (clk),
      .rst_n (rst_n),
      .valid_in (valid_out),
      .fft_imag_in (out_imag),
      .fft_imag_ref_in (fft_imag_ref_in),
      .symbol_index (write_addr),
      .comparison_result (imag_comparison_result)
    );

    FFT_out_imag u_imag_ref (
        .clka (clk),
        .addra (write_addr),
        .douta (fft_imag_ref_in)
    );

    FFT_output_real u_real_ref (
        .clka (clk),
        .addra (write_addr),
        .douta (fft_real_ref_in)
    );

    dds_counter #(.ADDR_WIDTH(ADDR_WIDTH + 1)) u_count (
        .enable (start),
        .clk (clk),
        .rst_n (rst_n),
        .count (dds_count)
    );


    
endmodule