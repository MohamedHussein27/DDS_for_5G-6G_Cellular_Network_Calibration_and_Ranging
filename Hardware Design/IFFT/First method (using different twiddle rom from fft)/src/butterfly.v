// Analog Devices 
// GP Ain-shams University
// Butterfly unit for 4096-FFT 


/* Description ..........
    butterfly unit is the basic building block of the FFT algorithm, it takes 2 complex inputs and produces 2 complex outputs,
    the first output is the sum of the 2 inputs and the second output is the difference of the 2 inputs, this is repeated 
    for log2(N) stages till the final result.
*/


module butterfly #(parameter WL = 16) (
    input wire signed [WL-1:0] in1_real, in1_imag, in2_real, in2_imag,
    output wire signed [WL-1:0] a_real, a_imag, b_real, b_imag
);
    assign a_real = in1_real + in2_real;  // a is the first "upper" output of the butterfly unit, which is the sum of the two inputs 
    assign a_imag = in1_imag + in2_imag;
    assign b_real = in1_real - in2_real;  // b is the second "lower" output of the butterfly unit, which is the difference of the two inputs
    assign b_imag = in1_imag - in2_imag;
endmodule