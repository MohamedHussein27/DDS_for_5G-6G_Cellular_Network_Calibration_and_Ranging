package ALU_driver_pkg;

  // Import ALU transaction definition
  import ALU_transaction_pkg::*;

  // ALU driver class
  // Responsible for driving stimulus onto the DUT via the interface
  class ALU_driver;

    // Mailbox for receiving transactions from generator
    mailbox #(ALU_transaction) gen_to_dr;

    // Transaction handle used by the driver
    ALU_transaction tr_dr;

    // Virtual interface handle (TB modport)
    virtual ALU_if.tb vif;

    // Constructor
    // - Receives virtual interface
    // - Receives mailbox from generator
    function new(
      virtual ALU_if.tb vif,
      mailbox #(ALU_transaction) gen_to_dr
    );
      this.vif        = vif;
      this.gen_to_dr  = gen_to_dr;
      tr_dr           = new();
    endfunction : new

    // Main driver task
    task run();
      forever begin

        // Get a transaction from the generator
        gen_to_dr.get(tr_dr);

        // Drive inputs onto the DUT
        vif.a  = tr_dr.a;
        vif.b  = tr_dr.b;
        vif.op = tr_dr.op;

        // Wait for one clock cycle (or synchronization point)
        @(negedge vif.clk);
      end
    endtask : run

  endclass : ALU_driver

endpackage : ALU_driver_pkg
