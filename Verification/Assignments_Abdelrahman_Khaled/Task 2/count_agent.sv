package count_agent_pkg;

  // Import required packages
  import count_transaction_pkg::*;
  import count_driver_pkg::*;
  import count_monitor_pkg::*;
  import count_generator_pkg::*;

  // ---------------------------------------------------------------------------
  // count_agent class
  // Encapsulates generator, driver, and monitor, along with mailboxes
  // ---------------------------------------------------------------------------
  class count_agent;

    // =========================
    // COMPONENTS
    // =========================
    count_driver    drv;   // Drives stimulus to DUT
    count_monitor   mon;   // Monitors DUT outputs
    count_generator gen;   // Generates transactions

    // =========================
    // VIRTUAL INTERFACES
    // =========================
    virtual count_if.tb      vif_tb;   // For driver
    virtual count_if.monitor vif_mon;  // For monitor

    // =========================
    // MAILBOXES
    // =========================
    mailbox #(count_transaction) mbx_gen_to_dr;   // Generator → Driver
    mailbox #(count_transaction) mbx_mon_to_sb;   // Monitor → Scoreboard
    mailbox #(count_transaction) mbx_mon_to_cov;  // Monitor → Coverage

    // =========================
    // CONSTRUCTOR
    // =========================
    function new();
    
      // Create components
      drv = new();
      mon = new();
      gen = new();

      // Create mailboxes
      mbx_gen_to_dr  = new();
      mbx_mon_to_sb  = new();
      mbx_mon_to_cov = new();
    endfunction : new

    // =========================
    // CONNECT FUNCTION
    // Assign interfaces and connect components to mailboxes
    // =========================
    function void connect(
      virtual count_if.tb      vif_tb,
      virtual count_if.monitor vif_mon
    );
      this.vif_tb  = vif_tb;
      this.vif_mon = vif_mon;

      // Instantiate components
      gen.connect(mbx_gen_to_dr);                          // Generator → Driver mailbox
      drv.connect(vif_tb, mbx_gen_to_dr);                  // Driver → DUT
      mon.connect(vif_mon, mbx_mon_to_sb, mbx_mon_to_cov); // Monitor → SB & coverage
    endfunction : connect

    // =========================
    // RUN PHASE
    // =========================
    task run();
      fork
        gen.run();
        drv.run();
        mon.run();
      join_none
    endtask : run

  endclass : count_agent

endpackage : count_agent_pkg
