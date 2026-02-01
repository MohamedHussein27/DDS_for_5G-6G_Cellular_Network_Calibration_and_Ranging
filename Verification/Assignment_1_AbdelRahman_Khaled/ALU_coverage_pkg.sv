package ALU_coverage_pkg;

  // Import ALU transaction definition
  import ALU_transaction_pkg::*;

  // Coverage class for ALU functional coverage
  class ALU_coverage;

    // Transaction handle used for coverage sampling
    ALU_transaction tr_cvg;

    // Covergroup to track ALU operation coverage
    covergroup ALU_cg;

      // Coverpoint for ALU operation field
      coverpoint tr_cvg.op {
        bins ADD = {0};
        bins XOR = {1};
        bins AND = {2};
        bins OR  = {3};
      }

      // Coverpoint for ALU operand a (4-bit)
      coverpoint tr_cvg.a {
          bins values_a[] = {[0:15]};  // 4-bit range 0..15
      }

      // Coverpoint for ALU operand b (4-bit)
      coverpoint tr_cvg.b {
          bins values_b[] = {[0:15]};  // 4-bit range 0..15
      }


    endgroup : ALU_cg

    // Constructor
    function new();
      tr_cvg = new();   // Create transaction object
      ALU_cg = new();   // Create covergroup
    endfunction : new

    // Sample function called by the monitor
    function void sample_data(ALU_transaction tr_DUT);
      // Copy transaction data to coverage object
      tr_cvg = tr_DUT;

      // Sample coverage
      ALU_cg.sample();
    endfunction : sample_data

  endclass : ALU_coverage

endpackage : ALU_coverage_pkg
