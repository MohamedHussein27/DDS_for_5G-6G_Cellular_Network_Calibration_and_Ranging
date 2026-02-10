package count_driver_pkg;

  // Import count transaction definition
  import count_transaction_pkg::*;

  // ------------------------------------------------------------
  // count Driver
  // Responsible for driving stimulus from the generator
  // onto the DUT through the count interface
  // ------------------------------------------------------------
  class count_driver;

    // Mailbox used to receive transactions from the generator
    mailbox #(count_transaction) gen_to_dr;

    // Transaction handle used by the driver
    count_transaction tr_dr;

    // Virtual interface handle (TB modport)
    // Used to drive signals into the DUT
    virtual count_if.tb vif;

    // ------------------------------------------------------------
    // Constructor
    // Kept intentionally empty
    // Object connections are handled explicitly in connect()
    // ------------------------------------------------------------
    function new();
      tr_dr = new();
    endfunction

    // ------------------------------------------------------------
    // Connect function
    // Binds the virtual interface and mailbox to the driver
    // ------------------------------------------------------------
    function void connect(
      virtual count_if.tb vif,
      mailbox #(count_transaction) gen_to_dr
    );
      this.vif       = vif;
      this.gen_to_dr = gen_to_dr;
    endfunction

    // ------------------------------------------------------------
    // Main driver task
    // Continuously receives transactions and drives them
    // onto the count DUT inputs
    // ------------------------------------------------------------
    task run();

      // Note:
      // Reset sequencing is assumed to be handled externally
      // (either by the test or the generator)

      forever begin

        // Block until a transaction is received
        gen_to_dr.get(tr_dr);

        // Drive transaction fields onto the interface signals
        vif.rst_n      = tr_dr.rst_n; 
        vif.start      = tr_dr.start;
        vif.flag       = tr_dr.flag;
        vif.wait_timer = tr_dr.wait_timer;

        // Synchronize to the clock
        // Inputs are held stable until the next cycle
        @(negedge vif.clk);

      end
    endtask

  endclass : count_driver

endpackage : count_driver_pkg
