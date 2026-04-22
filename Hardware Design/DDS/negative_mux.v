module negative_mux (
    input [7:0] amplitude_out,
    input negative_flag,
    output [7:0] final_amplitude
);
    
assign final_amplitude = (negative_flag) ? (-amplitude_out) : amplitude_out;

endmodule