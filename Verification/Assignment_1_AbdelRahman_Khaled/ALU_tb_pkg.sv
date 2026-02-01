package ALU_tb_pkg;

  // Import ALU transaction definition
  import ALU_transaction_pkg::*;

  // Import shared package (test_finished flag, sample event, etc.)
  import ALU_shared_pkg::*;

  // ALU driver / testbench class
  class ALU_tb;

    // Virtual interface handle (TB modport)
    virtual ALU_if.tb vif;

    // Transaction object used to generate randomized stimulus
    ALU_transaction tr;

    // Constructor
    function new(virtual ALU_if.tb vif);
      this.vif = vif;
      tr = new();
    endfunction : new

    // Main stimulus generation task
    task run();

      // Initialize DUT inputs
      vif.a  = '0;
      vif.b  = '0;
      vif.op = '0;
      test_finished = 0;

      // Apply a fixed number of randomized transactions
      repeat (1000) begin
        // Small delay between transactions
        #10;

        // Randomize ALU transaction fields (a, b, op)
        assert(tr.randomize());

        // Drive randomized values onto DUT inputs
        vif.a  = tr.a;
        vif.b  = tr.b;
        vif.op = tr.op;
      end

      // Indicate that the test has finished
      test_finished = 1;

    endtask : run

  endclass : ALU_tb

endpackage : ALU_tb_pkg
