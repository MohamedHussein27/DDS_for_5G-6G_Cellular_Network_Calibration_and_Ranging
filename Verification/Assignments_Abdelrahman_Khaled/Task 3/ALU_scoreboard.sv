package ALU_scoreboard_pkg;

  import uvm_pkg::*;
  `include "uvm_macros.svh"

  import ALU_seq_item_pkg::*;

  //----------------------------------------------------------------------------
  // ALU Scoreboard
  //----------------------------------------------------------------------------
  // Responsibilities:
  // 1) Maintain a reference model for ALU operations
  // 2) Compare DUT transactions to reference results
  // 3) Count errors and correct transactions
  // 4) Broadcast results for reporting/coverage
  //----------------------------------------------------------------------------
  class ALU_scoreboard extends uvm_scoreboard;

    `uvm_component_utils(ALU_scoreboard)

    // Reference signals
    logic       c_ref;
    logic [3:0] out_ref;

    // Operation enumeration
    typedef enum logic [1:0] {
      ADD_op = 0,
      XOR_op = 1,
      AND_op = 2,
      OR_op  = 3
    } operation_t;

    // Counters
    int Error_count   = 0;
    int Correct_count = 0;

    // Transaction handle
    ALU_seq_item item_sb;

    // TLM ports
    uvm_analysis_export #(ALU_seq_item) sb_export;
    uvm_tlm_analysis_fifo #(ALU_seq_item) sb_fifo;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_scoreboard", uvm_component parent = null);
      super.new(name, parent);
    endfunction

    //--------------------------------------------------------------------------
    // Build Phase
    //--------------------------------------------------------------------------
    function void build_phase(uvm_phase phase);
      super.build_phase(phase);
      sb_export = new("sb_export", this);
      sb_fifo   = new("sb_fifo", this);
    endfunction

    //--------------------------------------------------------------------------
    // Connect Phase
    //--------------------------------------------------------------------------
    function void connect_phase(uvm_phase phase);
      super.connect_phase(phase);
      sb_export.connect(sb_fifo.analysis_export);
    endfunction

    //--------------------------------------------------------------------------
    // Reference Model
    //--------------------------------------------------------------------------
    function void ALU_ref(ALU_seq_item item);
      unique case(item.op)
        ADD_op: {c_ref, out_ref} = item.a + item.b;
        XOR_op: begin out_ref = item.a ^ item.b; c_ref = 0; end
        AND_op: begin out_ref = item.a & item.b; c_ref = 0; end
        OR_op : begin out_ref = item.a | item.b; c_ref = 0; end
        default: begin out_ref = '0; c_ref = 0; end
      endcase
    endfunction

    //--------------------------------------------------------------------------
    // Compare DUT transaction with reference
    //--------------------------------------------------------------------------
    function void ALU_check(ALU_seq_item item);
      ALU_ref(item);

      if ((item.op == ADD_op && {item.c, item.out} !== {c_ref, out_ref}) ||
          (item.op != ADD_op && (item.out !== out_ref || item.c !== 0))) begin
        Error_count++;
        `uvm_error("ALU_CHECK", 
          $sformatf("Transaction mismatch! DUT: %s, REF: {c=%0b, out=%0b}", 
                    item.convert2string(), c_ref, out_ref))
      end else begin
        Correct_count++;
      end
    endfunction

    //--------------------------------------------------------------------------
    // Run Phase
    //--------------------------------------------------------------------------
    task run_phase(uvm_phase phase);
      super.run_phase(phase);

      forever begin
        // Get transactions from FIFO (connected to monitor)
        sb_fifo.get(item_sb);

        // Compare with reference model
        ALU_check(item_sb);
      end
    endtask

    //--------------------------------------------------------------------------
    // Report Phase
    //--------------------------------------------------------------------------
    function void report_phase(uvm_phase phase);
      super.report_phase(phase);
      `uvm_info("ALU_SCOREBOARD", 
                $sformatf("Total correct transactions: %0d", Correct_count), UVM_MEDIUM)
      `uvm_info("ALU_SCOREBOARD", 
                $sformatf("Total error transactions:   %0d", Error_count), UVM_MEDIUM)
    endfunction

  endclass : ALU_scoreboard

endpackage : ALU_scoreboard_pkg
