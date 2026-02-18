package alu_scoreboard_pkg;

  import uvm_pkg::*;
  import alu_pkg::*;
  import alu_sequence_item_pkg::*;
  `include "uvm_macros.svh"

  class alu_scoreboard extends uvm_scoreboard;

    `uvm_component_utils(alu_scoreboard)

    // --------------------------------------------
    // 1. Define Struct for Error Tracking
    // --------------------------------------------
    
    // --------------------------------------------
    // 2. Queue to store failed transactions
    // --------------------------------------------
    error_packet_t error_log_q[$];

    // --------------------------------------------
    // TLM Connections
    // --------------------------------------------
    uvm_analysis_export #(alu_sequence_item) sb_export;
    uvm_tlm_analysis_fifo #(alu_sequence_item) sb_fifo;

    // --------------------------------------------
    // Counters
    // --------------------------------------------
    int error_count   = 0;
    int correct_count = 0;

    // --------------------------------------------
    // Constructor
    // --------------------------------------------
    function new(string name = "alu_scoreboard",
                 uvm_component parent = null);
      super.new(name, parent);
    endfunction


    // --------------------------------------------
    // Build Phase
    // --------------------------------------------
    function void build_phase(uvm_phase phase);
      super.build_phase(phase);

      sb_export = new("sb_export", this);
      sb_fifo   = new("sb_fifo", this);

    endfunction


    // --------------------------------------------
    // Connect Phase
    // --------------------------------------------
    function void connect_phase(uvm_phase phase);
      super.connect_phase(phase);

      sb_export.connect(sb_fifo.analysis_export);

    endfunction


    // --------------------------------------------
    // Run Phase
    // --------------------------------------------
    task run_phase(uvm_phase phase);

      
      alu_sequence_item item;
      

      super.run_phase(phase);

      forever begin

        sb_fifo.get(item);

        // ----------------------------
        // RESET HANDLING
        // ----------------------------
        if (!item.rst_n) begin
          expected = 5'b0;
        end
        else begin
           case (item.op)
            ADD_op: expected = item.a + item.b;
            XOR_op: expected = item.a ^ item.b;
            AND_op: expected = item.a & item.b;
            OR_op : expected = item.a | item.b;
            default: expected = 5'b0;
          endcase
        end

        // ----------------------------
        // COMPARE
        // ----------------------------
        if (item.result !== expected) begin
          error_count++;

          // 3. CAPTURE THE ERROR
          err_pkt.a  = item.a;
          err_pkt.b  = item.b;
          err_pkt.op = item.op;
          error_log_q.push_back(err_pkt); // Save to queue

          `uvm_error("ALU_SCOREBOARD",
            $sformatf("ERROR: A=%0d B=%0d OP=%s EXP=%0d GOT=%0d | errors=%0d",
              item.a, item.b, item.op.name(), // Using .name() for enum string
              expected, item.result, error_count))

        end
        else begin
          correct_count++;

          `uvm_info("ALU_SCOREBOARD",
            $sformatf("CORRECT: A=%0d B=%0d OP=%s EXP=%0d GOT=%0d | correct=%0d",
              item.a, item.b, item.op.name(),
              expected, item.result, correct_count),
            UVM_LOW)
        end

      end

    endtask


    // --------------------------------------------
    // Report Phase
    // --------------------------------------------
    function void report_phase(uvm_phase phase);
      error_packet_t unique_errors_q[$]; // Local queue for filtering

      super.report_phase(phase);

      `uvm_info("ALU_REPORT",
        $sformatf("Total successful transactions: %0d",
                  correct_count),
        UVM_LOW)

      `uvm_info("ALU_REPORT",
        $sformatf("Total failed transactions: %0d",
                  error_count),
        UVM_LOW)

      // 4. PRINT UNIQUE ERRORS
      if (error_count > 0) begin
        // Filter duplicates
        unique_errors_q = error_log_q.unique(); 

        `uvm_info("ALU_REPORT", "========================================", UVM_LOW)
        `uvm_info("ALU_REPORT", $sformatf(" UNIQUE FAILURE PATTERNS FOUND: %0d", unique_errors_q.size()), UVM_LOW)
        `uvm_info("ALU_REPORT", "========================================", UVM_LOW)

        foreach(unique_errors_q[i]) begin
           `uvm_info("ALU_REPORT", 
              $sformatf("Pattern #%0d:  Op=%s  A=%0d  B=%0d", 
              i+1, unique_errors_q[i].op.name(), unique_errors_q[i].a, unique_errors_q[i].b), 
              UVM_LOW)
        end
        `uvm_info("ALU_REPORT", "========================================", UVM_LOW)
        
        `uvm_error("ALU_REPORT", "TEST FAILED ")
      end
      else begin
        `uvm_info("ALU_REPORT", "TEST PASSED ", UVM_NONE)
      end

    endfunction

  endclass

endpackage