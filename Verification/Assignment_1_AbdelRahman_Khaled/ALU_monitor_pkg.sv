package ALU_monitor_pkg;

  // Import required packages
  import ALU_transaction_pkg::*;
  import ALU_scoreboard_pkg::*;
  import ALU_coverage_pkg::*;
  import ALU_shared_pkg::*;

  // ALU Monitor class
  class ALU_monitor;

    // Virtual interface handle (monitor modport)
    virtual ALU_if.monitor vif;

    // Transaction object to hold sampled DUT values
    ALU_transaction tr;

    // Scoreboard instance
    ALU_scoreboard tr_sb;

    // Coverage instance
    ALU_coverage tr_cov;

    // Constructor
    function new(virtual ALU_if.monitor vif);
      this.vif = vif;
      tr     = new();
      tr_sb  = new();
      tr_cov = new();
    endfunction : new

    // Main monitoring task
    task run();

      #5;
      forever begin

        #10; // Wait for stimulus to stabilize

        // Sample ALU inputs and outputs from the interface
        tr.a   = vif.a;
        tr.b   = vif.b;
        tr.op  = vif.op;
        tr.out = vif.out;
        tr.c   = vif.c;

        // Send the sampled transaction to the scoreboard
        tr_sb.ALU_check(tr);

        // Send the sampled transaction to the coverage collector
        tr_cov.sample_data(tr);

        // Check for end-of-test condition
        if (test_finished == 1) begin
          $display("Error_count=%0d, Correct_count=%0d",
                   Error_count, Correct_count);
          $stop;
        end
      end
    endtask : run

  endclass : ALU_monitor

endpackage : ALU_monitor_pkg
