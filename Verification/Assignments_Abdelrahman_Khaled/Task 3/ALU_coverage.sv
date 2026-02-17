package ALU_coverage_pkg;

  import uvm_pkg::*;
  `include "uvm_macros.svh"

  import ALU_seq_item_pkg::*;

  //----------------------------------------------------------------------------
  // ALU Coverage
  //----------------------------------------------------------------------------
  // Responsibilities:
  // 1) Collect transactions from monitor
  // 2) Sample data into covergroups
  // 3) Provide analysis port for connection to monitor
  //----------------------------------------------------------------------------
  class ALU_coverage extends uvm_component;

    `uvm_component_utils(ALU_coverage)

    // Transaction handle
    ALU_seq_item item_cov;

    // TLM ports
    uvm_analysis_export #(ALU_seq_item) cov_export;
    uvm_tlm_analysis_fifo #(ALU_seq_item) cov_fifo;

    //--------------------------------------------------------------------------
    // ALU Functional Covergroup
    //--------------------------------------------------------------------------
    covergroup ALU_cg;

      // Coverpoint for ALU operation
      coverpoint item_cov.op {
        bins ADD = {0};
        bins XOR = {1};
        bins AND = {2};
        bins OR  = {3};
      }

      // Coverpoint for operand a (4-bit)
      coverpoint item_cov.a {
        bins values_a[] = {[0:15]};
      }

      // Coverpoint for operand b (4-bit)
      coverpoint item_cov.b {
        bins values_b[] = {[0:15]};
      }

    endgroup : ALU_cg

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_coverage", uvm_component parent = null);
      super.new(name, parent);

      // Instantiate covergroup
      ALU_cg = new();
    endfunction

    //--------------------------------------------------------------------------
    // Build Phase
    //--------------------------------------------------------------------------
    function void build_phase(uvm_phase phase);
      super.build_phase(phase);

      // Create TLM ports
      cov_export = new("cov_export", this);
      cov_fifo   = new("cov_fifo", this);
    endfunction

    //--------------------------------------------------------------------------
    // Connect Phase
    //--------------------------------------------------------------------------
    function void connect_phase(uvm_phase phase);
      super.connect_phase(phase);
      cov_export.connect(cov_fifo.analysis_export);
    endfunction

    //--------------------------------------------------------------------------
    // Sample data into covergroup
    //--------------------------------------------------------------------------
    function void sample_data(ALU_seq_item item);
      item_cov = item;
      ALU_cg.sample();
    endfunction

    //--------------------------------------------------------------------------
    // Run phase: continuously get transactions from FIFO and sample
    //--------------------------------------------------------------------------
    task run_phase(uvm_phase phase);
      super.run_phase(phase);

      forever begin
        // Get transaction from FIFO (connected to monitor)
        cov_fifo.get(item_cov);

        // Sample transaction into covergroup
        sample_data(item_cov);
      end
    endtask

  endclass : ALU_coverage

endpackage : ALU_coverage_pkg
