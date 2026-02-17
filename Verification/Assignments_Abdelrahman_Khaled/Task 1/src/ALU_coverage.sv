package ALU_coverage_pkg;

  import ALU_transaction_pkg::*;

  // ALU functional coverage class
  class ALU_coverage;

    // Transaction handle for sampling
    ALU_transaction tr_cvg;

    // Mailbox connected to monitor
    mailbox #(ALU_transaction) mbx_mon_to_cov;

    // Covergroup for ALU functional coverage
    covergroup ALU_cg;

      // Coverpoint for ALU operation
      coverpoint tr_cvg.op {
        bins ADD = {0};
        bins XOR = {1};
        bins AND = {2};
        bins OR  = {3};
      }

      // Coverpoint for operand a (4-bit)
      coverpoint tr_cvg.a {
        bins values_a[] = {[0:15]};
      }

      // Coverpoint for operand b (4-bit)
      coverpoint tr_cvg.b {
        bins values_b[] = {[0:15]};
      }

    endgroup : ALU_cg

    // Constructor: initializes mailbox, transaction object, and covergroup
    function new(mailbox #(ALU_transaction) mbx_mon_to_cov);
      this.mbx_mon_to_cov = mbx_mon_to_cov;
      tr_cvg = new();       // Transaction object for covergroup
      ALU_cg = new();       // Create the covergroup
    endfunction : new

    // Sample a transaction into the covergroup
    function void sample_data(ALU_transaction tr_DUT);
      tr_cvg.a  = tr_DUT.a;
      tr_cvg.b  = tr_DUT.b;
      tr_cvg.op = tr_DUT.op;

      ALU_cg.sample();
    endfunction : sample_data

    // Main run task: continuously sample transactions from monitor
    task run();
      forever begin
        mbx_mon_to_cov.get(tr_cvg);   // Get new transaction from monitor
        sample_data(tr_cvg);          // Sample it into the covergroup
      end
    endtask

  endclass : ALU_coverage

endpackage : ALU_coverage_pkg
