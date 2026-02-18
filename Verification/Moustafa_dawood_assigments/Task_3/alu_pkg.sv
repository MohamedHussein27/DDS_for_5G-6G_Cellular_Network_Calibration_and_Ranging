package alu_pkg;



typedef enum bit [1:0] {
    ADD_op=0,
     XOR_op=1,
      AND_op=2,
       OR_op=3
  } opcode_t;


 typedef struct packed {
      bit [3:0] a;
      bit [3:0] b;
      opcode_t  op;
    } error_packet_t;


      
      bit [4:0] expected;
      error_packet_t err_pkt; // Temp struct variable

      

    
endpackage