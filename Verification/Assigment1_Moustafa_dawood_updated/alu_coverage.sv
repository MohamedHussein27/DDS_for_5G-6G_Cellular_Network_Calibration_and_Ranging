
/*
class coverage ;

virtual alu_bfm bfm;

        bit [3:0] a;
        bit [3:0] b;
        opcode_t  op_set;

covergroup op_cov;

    option.per_instance = 1;

    // Coverpoint for operand A
    cp_operand_a : coverpoint a {
        bins low_values = {0, 1, 2, 3, 4, 5, 6, 7};
        bins high_values = {8, 9, 10, 11, 12, 13, 14, 15};
    }

    // Coverpoint for operand B
    cp_operand_b : coverpoint b {
        bins low_values = {0, 1, 2, 3, 4, 5, 6, 7};
        bins high_values = {8, 9, 10, 11, 12, 13, 14, 15};
    }

    // Coverpoint for operation code
    cp_opcode : coverpoint op_set {
        bins add_op = {ADD_op};
        bins xor_op = {XOR_op};
        bins and_op = {AND_op};
        bins or_op = {OR_op};
    }




endgroup : op_cov


covergroup alu_corner_cases_cg;
    // Coverpoint for operations
    all_ops : coverpoint op_set {
        bins add_op = {0};  // ADD_op = 0
        bins xor_op = {1};  // XOR_op = 1  
        bins and_op = {2};  // AND_op = 2
        bins or_op  = {3};  // OR_op = 3
        illegal_bins invalid = default; // Should never get other values
    }

    // Coverpoint for input A (4-bit)
    a_input : coverpoint a {
        bins zero    = {4'b0000};           // All zeros
        bins one     = {4'b0001};           // 1 (min+1)
        bins max_minus_one = {4'b1110};     // 14 
        bins max     = {4'b1111};           // All ones
        bins others  = {[4'b0010:4'b1101]}; // All other values 2-13
    }

    // Coverpoint for input B (4-bit)  
    b_input : coverpoint b {
        bins zero    = {4'b0000};           // All zeros
        bins one     = {4'b0001};           // 1 
        bins max_minus_one = {4'b1110};     // 14 
        bins max     = {4'b1111};           // All ones
        bins others  = {[4'b0010:4'b1101]}; // All other values 2-13
    }

    // CROSS COVERAGE 
    corner_cross : cross all_ops, a_input, b_input {
        // Test ZEROS with each operation
        bins add_zero_zero = binsof(all_ops.add_op) && 
                            (binsof(a_input.zero) && binsof(b_input.zero));
        
        bins add_zero_any = binsof(all_ops.add_op) &&
                           (binsof(a_input.zero) || binsof(b_input.zero));
        
        bins xor_zero_zero = binsof(all_ops.xor_op) && 
                            (binsof(a_input.zero) && binsof(b_input.zero));
        
        bins and_zero_zero = binsof(all_ops.and_op) && 
                            (binsof(a_input.zero) && binsof(b_input.zero));
        
        bins or_zero_zero = binsof(all_ops.or_op) && 
                           (binsof(a_input.zero) && binsof(b_input.zero));
        
        // Test ONES (max) with each operation  
        bins add_max_max = binsof(all_ops.add_op) &&
                          (binsof(a_input.max) && binsof(b_input.max));
        
        bins add_max_any = binsof(all_ops.add_op) &&
                          (binsof(a_input.max) || binsof(b_input.max));
        
        bins xor_max_max = binsof(all_ops.xor_op) &&
                          (binsof(a_input.max) && binsof(b_input.max));
        
        bins and_max_max = binsof(all_ops.and_op) &&
                          (binsof(a_input.max) && binsof(b_input.max));
        
        bins or_max_max = binsof(all_ops.or_op) &&
                         (binsof(a_input.max) && binsof(b_input.max));
        
        // Test critical edge: max-1 values (14)
        bins add_near_overflow = binsof(all_ops.add_op) &&
                                ((binsof(a_input.max_minus_one) && binsof(b_input.one)) ||
                                 (binsof(a_input.one) && binsof(b_input.max_minus_one)) ||
                                 (binsof(a_input.max_minus_one) && binsof(b_input.max_minus_one)));
        
        // Test operation transitions (same inputs, different all_ops)
        bins same_input_diff_ops = 
            (binsof(a_input.zero) && binsof(b_input.zero)) ||
            (binsof(a_input.max) && binsof(b_input.max)) ||
            (binsof(a_input.one) && binsof(b_input.one));
        
        
        // Specifically track carry generation for ADD
        bins add_with_carry = binsof(all_ops.add_op) && 
                             ((a + b) > 15);  // Will generate carry
    }
endgroup


function new (virtual alu_bfm b);
     op_cov = new();
     alu_corner_cases_cg = new();
     bfm = b;
   endfunction : new



    task sample();
    forever begin
        @(negedge bfm.clk)

          a = bfm.operand_a;
          b = bfm.operand_b;
          op_set  = bfm.op_set;
          op_cov.sample();
          alu_corner_cases_cg.sample();
    end
          
    endtask : sample
endclass

*/
 /*
class alu_coverage;
 

  virtual alu_bfm bfm;

  bit [3:0] a, b;
  opcode_t op;

  covergroup alu_cg;
    cp_a  : coverpoint a;
    cp_b  : coverpoint b;
    cp_op : coverpoint op;
    cross cp_a, cp_b, cp_op;
  endgroup

  function new(virtual alu_bfm b);
    bfm = b;
    alu_cg = new();
  endfunction

   task run();
    forever begin
      @(posedge bfm.clk);  // sample at positive edge
      a  = bfm.operand_a;
      b  = bfm.operand_b;
      op = bfm.op_set;
      alu_cg.sample();
    end
  endtask

endclass
*/
//import alu_pkg::*;
package alu_coverage_pkg;
import alu_pkg::*;

    

class alu_coverage;

  virtual alu_bfm bfm;
  bit [3:0] a, b;
  opcode_t op;

  covergroup alu_cg;
    cp_a  : coverpoint a;
    cp_b  : coverpoint b;
    cp_op : coverpoint op;
    cross cp_a, cp_b, cp_op;
  endgroup

  function new();
  alu_cg = new();
  endfunction

  task run();
    forever begin
      @(posedge bfm.clk);
      a  = bfm.operand_a;
      b  = bfm.operand_b;
      op = bfm.op_set;
      alu_cg.sample();
    end
  endtask

endclass
endpackage