package alu_coverage_pkg;
  import uvm_pkg::*;
  import alu_pkg::*;
  import alu_sequence_item_pkg::*; 
  `include "uvm_macros.svh"

  // 1. Extend uvm_component (Since we are managing the FIFO ourselves)
  class alu_coverage extends uvm_component;
    
    `uvm_component_utils(alu_coverage)

    //  The Components
    // The "Mailbox" (Internal Buffer)
    uvm_tlm_analysis_fifo #(alu_sequence_item) cov_fifo; 
    
    // The  (External Port) - Monitor connects here
    uvm_analysis_export #(alu_sequence_item) analysis_export; 

    // Local variables for coverage sampling
    bit [3:0] a, b;
    opcode_t  op;
    alu_sequence_item item; // Temp variable for reading

    // 3. Covergroup
    covergroup alu_cg;
     // option.per_instance = 1; 
      cp_a  : coverpoint a;
      cp_b  : coverpoint b;
      cp_op : coverpoint op;
      cross cp_a, cp_b, cp_op;
    endgroup
    
    // 4. Constructor
    function new(string name, uvm_component parent);
      super.new(name, parent);
      // Create the covergroup
      alu_cg = new(); 
      
      // Create the FIFO and Export
      cov_fifo = new("cov_fifo", this);
      analysis_export = new("analysis_export", this);
    endfunction

    // 5. Connect Phase
    // This connects the outside world (Monitor) to our internal FIFO
    function void connect_phase(uvm_phase phase);
      super.connect_phase(phase);
      analysis_export.connect(cov_fifo.analysis_export);
    endfunction

    // 6. Run Phase
    
    task run_phase(uvm_phase phase);
      super.run_phase(phase);
      
      forever begin
        // BLOCKING: Wait for data to arrive in FIFO
        cov_fifo.get(item); 
        
        // Copy data to local variables
        a  = item.a;
        b  = item.b;
        op = item.op;
        
        // Sample
        alu_cg.sample();
      end
    endtask

  endclass
endpackage