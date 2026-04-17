// Analog Devices 
// GP Ain-shams University
// Delay feedback unit for 4096-FFT 


/* Description ..........
   L = 2^(n-s) where n is the log 2 (FFT size = 4096) and s is the stage number.
   we used a shift register in this block to store the first half of the inputs in the register 
   and then when the 2nd half come we shift out the stored value one by one to pass it with 
   the incoming input to the butterfly to do it's operation, and the delay depend on the stage number 
   for Ex: 
   stage 12 we need 2048 delay and for stage 11 we need 1024 and so on ...
*/



module delayfeedback #(
    parameter L = 2048, // Delay length (depends on the stage number, for stage 12 it's 2048)
    parameter WL = 16    // Word length
)(
    input wire clk,
    input wire rst_n,
    input wire signed [WL-1:0] data_in_real,
    input wire signed [WL-1:0] data_in_imag,
    output wire signed [WL-1:0] data_out_real,
    output wire signed [WL-1:0] data_out_imag
);


reg [WL-1:0] delay_reg_real [L-1:0]; // Shift register for delay real part
reg [WL-1:0] delay_reg_imag [L-1:0]; // Shift register for delay imaginary part   
integer i=0,j=0; // Loop variables   

always @(posedge clk, negedge rst_n) begin 
    if(!rst_n) begin
        for (i= 0; i < L; i = i + 1) begin
            delay_reg_real[i] <= 0;
            delay_reg_imag[i] <= 0;
        end
    end else begin
        for (j = L-1; j > 0; j = j - 1) begin
            delay_reg_real[j] <= delay_reg_real[j-1]; // Shift the register  >> Example [3 5 7 8 1] [5 7 8 1 data_in] dataout=3.
            delay_reg_imag[j] <= delay_reg_imag[j-1];
        end
        delay_reg_real[0] <= data_in_real; // Load new input at the beginning of the register
        delay_reg_imag[0] <= data_in_imag;
    end
    end
assign data_out_real = delay_reg_real[L-1];    
assign data_out_imag = delay_reg_imag[L-1];
endmodule