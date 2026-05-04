module ifft_4096_top #(
    parameter WL = 16
)(
    input wire clk,
    input wire rst_n,
    input wire valid_in,
    input wire signed [WL-1:0] in_real,
    input wire signed [WL-1:0] in_imag,
    output wire valid_out,
    output reg signed [WL-1:0] out_real,
    output reg signed [WL-1:0] out_imag
);
    wire [11:0] global_sel;
    wire [143:0] global_addr_bus;
    wire pipeline_en;

    controlunit_4096 #(.N(4096)) ctrl (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(valid_in),
        .sel(global_sel),
        .addr_bus(global_addr_bus),
        .pipeline_en(pipeline_en),
        .valid_out(valid_out)
    );

    wire signed [WL-1:0] stage_re [0:12];
    wire signed [WL-1:0] stage_im [0:12];

    assign stage_re[0] = in_real;
    assign stage_im[0] = in_imag;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            out_real <= 0;
            out_imag <= 0;
        end else begin
            out_real <= stage_re[12];
            out_imag <= stage_im[12];
        end
    end

    // --- STAGE 1 ---
    wire [(2*WL-1):0] twiddle_value_stg1;
    Twiddle_stage_1 rom_s1 (.clka(clk), .addra(global_addr_bus[10:0]), .douta(twiddle_value_stg1));

//   Twiddle_stage_1 rom_s1 (.clk(clk), .a(global_addr_bus[10:0]), .qspo(twiddle_value_stg1));

    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(2048)) stg1 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[11]), 
        .twiddle_re($signed(twiddle_value_stg1[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg1[15:0])),
        .in_real(stage_re[0]), .in_imag(stage_im[0]), 
        .out_real(stage_re[1]), .out_imag(stage_im[1])
    );

    // --- STAGE 2 ---
    wire [(2*WL-1):0] twiddle_value_stg2;
    Twiddle_Stage_2 rom_s2 (.clka(clk), .addra(global_addr_bus[21:12]), .douta(twiddle_value_stg2));
//    Twiddle_Stage_2 rom_s2 (.clk(clk), .a(global_addr_bus[21:12]), .qspo(twiddle_value_stg2));

    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(1024)) stg2 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[10]), 
        .twiddle_re($signed(twiddle_value_stg2[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg2[15:0])),
        .in_real(stage_re[1]), .in_imag(stage_im[1]), 
        .out_real(stage_re[2]), .out_imag(stage_im[2])
    );

    // --- STAGE 3 ---
    wire [(2*WL-1):0] twiddle_value_stg3;
    Twiddle_Stage_3 rom_s3 (.clka(clk), .addra(global_addr_bus[32:24]), .douta(twiddle_value_stg3));
    
    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(512)) stg3 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[9]), 
        .twiddle_re($signed(twiddle_value_stg3[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg3[15:0])),
        .in_real(stage_re[2]), .in_imag(stage_im[2]), 
        .out_real(stage_re[3]), .out_imag(stage_im[3])
    );

    // --- STAGE 4 ---
    wire [(2*WL-1):0] twiddle_value_stg4;
    Twiddle_Stage_4 rom_s4 (.clka(clk), .addra(global_addr_bus[43:36]), .douta(twiddle_value_stg4));

    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(256)) stg4 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[8]), 
        .twiddle_re($signed(twiddle_value_stg4[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg4[15:0])),
        .in_real(stage_re[3]), .in_imag(stage_im[3]), 
        .out_real(stage_re[4]), .out_imag(stage_im[4])
    );

    // --- STAGE 5 ---
    wire [(2*WL-1):0] twiddle_value_stg5;
//    Twiddle_Stage_5 rom_s5 (.clka(clk), .addra(global_addr_bus[54:48]), .douta(twiddle_value_stg5));
    Twiddle_Stage_5 rom_s5 (.clk(clk), .a(global_addr_bus[54:48]), .qspo(twiddle_value_stg5));
    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(128)) stg5 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[7]), 
        .twiddle_re($signed(twiddle_value_stg5[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg5[15:0])),
        .in_real(stage_re[4]), .in_imag(stage_im[4]), 
        .out_real(stage_re[5]), .out_imag(stage_im[5])
    );

    // --- STAGE 6 ---
    wire [(2*WL-1):0] twiddle_value_stg6;
//    Twiddle_Stage_6 rom_s6 (.clka(clk), .addra(global_addr_bus[65:60]), .douta(twiddle_value_stg6));
    Twiddle_Stage_6 rom_s6 (.clk(clk), .a(global_addr_bus[65:60]), .qspo(twiddle_value_stg6));

    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(64)) stg6 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[6]), 
        .twiddle_re($signed(twiddle_value_stg6[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg6[15:0])),
        .in_real(stage_re[5]), .in_imag(stage_im[5]), 
        .out_real(stage_re[6]), .out_imag(stage_im[6])
    );

    // --- STAGE 7 ---
    wire [(2*WL-1):0] twiddle_value_stg7;
   // Twiddle_Stage_7 rom_s7 (.clka(clk), .addra(global_addr_bus[76:72]), .douta(twiddle_value_stg7));
Twiddle_Stage_7 rom_s7 (.clk(clk), .a(global_addr_bus[76:72]), .qspo(twiddle_value_stg7));
    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(32)) stg7 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[5]), 
        .twiddle_re($signed(twiddle_value_stg7[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg7[15:0])),
        .in_real(stage_re[6]), .in_imag(stage_im[6]), 
        .out_real(stage_re[7]), .out_imag(stage_im[7])
    );

    // --- STAGE 8 ---
    wire [(2*WL-1):0] twiddle_value_stg8;
//    Twiddle_Stage_8 rom_s8 (.clka(clk), .addra(global_addr_bus[87:84]), .douta(twiddle_value_stg8));
Twiddle_Stage_8 rom_s8 (.clk(clk), .a(global_addr_bus[87:84]), .qspo(twiddle_value_stg8));
    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(16)) stg8 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[4]), 
        .twiddle_re($signed(twiddle_value_stg8[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg8[15:0])),
        .in_real(stage_re[7]), .in_imag(stage_im[7]), 
        .out_real(stage_re[8]), .out_imag(stage_im[8])
    );

    // --- STAGE 9 ---
    wire [(2*WL-1):0] twiddle_value_stg9;
//    Twiddle_Stage_9 rom_s9 (.clka(clk), .addra(global_addr_bus[98:96]), .douta(twiddle_value_stg9));
Twiddle_Stage_9 rom_s9 (.clk(clk), .a({1'b0,global_addr_bus[98:96]}), .qspo(twiddle_value_stg9));
    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(8)) stg9 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[3]), 
        .twiddle_re($signed(twiddle_value_stg9[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg9[15:0])),
        .in_real(stage_re[8]), .in_imag(stage_im[8]), 
        .out_real(stage_re[9]), .out_imag(stage_im[9])
    );

    // --- STAGE 10 ---
    wire [(2*WL-1):0] twiddle_value_stg10;
//    Twiddle_Stage_10 rom_s10 (.clka(clk), .addra(global_addr_bus[109:108]), .douta(twiddle_value_stg10));
Twiddle_Stage_10 rom_s10 (.clk(clk), .a({1'b0,global_addr_bus[109:108]}), .qspo(twiddle_value_stg10));
    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(4)) stg10 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[2]), 
        .twiddle_re($signed(twiddle_value_stg10[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg10[15:0])),
        .in_real(stage_re[9]), .in_imag(stage_im[9]), 
        .out_real(stage_re[10]), .out_imag(stage_im[10])
    );

    // --- STAGE 11 ---
    wire [(2*WL-1):0] twiddle_value_stg11;
 //   Twiddle_Stage_11 rom_s11 (.clka(clk), .addra(global_addr_bus[120]), .douta(twiddle_value_stg11));
Twiddle_Stage_11 rom_s11 (.clk(clk), .a({2'b00,global_addr_bus[120]}), .qspo(twiddle_value_stg11));
    ifft_stage_4096 #(.WL(WL), .DELAY_LEN(2)) stg11 (
        .clk(clk), .rst_n(rst_n), .valid_in(pipeline_en),
        .sel(global_sel[1]), 
        .twiddle_re($signed(twiddle_value_stg11[(2*WL-1):16])), 
        .twiddle_im($signed(twiddle_value_stg11[15:0])),
        .in_real(stage_re[10]), .in_imag(stage_im[10]), 
        .out_real(stage_re[11]), .out_imag(stage_im[11])
    );

    // --- STAGE 12 ---
    ifft_stage_last #(.WL(WL)) stg12 (
        .clk(clk), 
        .rst_n(rst_n), 
        .valid_in(pipeline_en),
        .sel(global_sel[0]), 
        .in_real(stage_re[11]), 
        .in_imag(stage_im[11]), 
        .out_real(stage_re[12]), 
        .out_imag(stage_im[12])
    );
endmodule