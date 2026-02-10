package ALU_agent_pkg;

  // Import required packages
  import ALU_transaction_pkg::*;
  import ALU_driver_pkg::*;
  import ALU_monitor_pkg::*;
  import ALU_generator_pkg::*;

  // ALU Agent class
  // Encapsulates generator, driver, and monitor
  class ALU_agent;

    // =========================
    // COMPONENTS
    // =========================
    ALU_driver    drv;   // Drives stimulus to DUT
    ALU_monitor   mon;   // Monitors DUT outputs
    ALU_generator gen;   // Generates transactions

    // =========================
    // VIRTUAL INTERFACES
    // =========================
    virtual ALU_if.tb      vif_tb;   // For driver
    virtual ALU_if.monitor vif_mon;  // For monitor

    // =========================
    // MAILBOXES
    // =========================
    mailbox #(ALU_transaction) mbx_gen_to_dr;   // Generator → Driver
    mailbox #(ALU_transaction) mbx_mon_to_sb;   // Monitor → Scoreboard
    mailbox #(ALU_transaction) mbx_mon_to_cov;  // Monitor → Coverage

    // =========================
    // CONSTRUCTOR
    // =========================
    function new(
      virtual ALU_if.tb      vif_tb,
      virtual ALU_if.monitor vif_mon
    );
      this.vif_tb  = vif_tb;
      this.vif_mon = vif_mon;

      // Create mailboxes inside the agent
      mbx_gen_to_dr  = new();
      mbx_mon_to_sb  = new();
      mbx_mon_to_cov = new();
    endfunction : new

    // =========================
    // BUILD PHASE
    // =========================
    function void build();

      // Generator: produces transactions
      gen = new(mbx_gen_to_dr);

      // Driver: consumes transactions and drives DUT
      drv = new(vif_tb, mbx_gen_to_dr);

      // Monitor: observes DUT and sends data to SB & coverage
      mon = new(vif_mon, mbx_mon_to_sb, mbx_mon_to_cov);

    endfunction : build

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

  endclass : ALU_agent

endpackage : ALU_agent_pkg
