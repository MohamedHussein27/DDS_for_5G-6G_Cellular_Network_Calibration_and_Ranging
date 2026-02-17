package ALU_env_pkg;

  // Import ALU transaction definition
  import ALU_transaction_pkg::*;
  import ALU_agent_pkg::*;
  import ALU_coverage_pkg::*;
  import ALU_scoreboard_pkg::*;

  class ALU_env;

    // =========================
    // COMPONENTS
    // =========================
    ALU_agent       ag;   // Agent encapsulating gen/driver/monitor
    ALU_scoreboard  sb;   // Scoreboard
    ALU_coverage    cov;  // Coverage

    // =========================
    // MAILBOXES
    // =========================
    mailbox #(ALU_transaction) mbx_mon_to_sb;
    mailbox #(ALU_transaction) mbx_mon_to_cov;

    // =========================
    // VIRTUAL INTERFACES
    // =========================
    virtual ALU_if.tb      vif_tb;   // Driver
    virtual ALU_if.monitor vif_mon;  // Monitor

    // =========================
    // CONSTRUCTOR
    // =========================
    function new(
      virtual ALU_if.tb      vif_tb,
      virtual ALU_if.monitor vif_mon
    );
      this.vif_tb  = vif_tb;
      this.vif_mon = vif_mon;

      // Create mailboxes
      mbx_mon_to_sb  = new();
      mbx_mon_to_cov = new();
    endfunction : new

    // =========================
    // BUILD PHASE
    // =========================
    function void build();

      // Create agent (it will create its internal mailbox: gen->drv)
      ag = new(vif_tb, vif_mon);
      ag.build();

      // Create scoreboard and coverage
      sb  = new(mbx_mon_to_sb);
      cov = new(mbx_mon_to_cov);

      // Connect agent monitor to scoreboard & coverage
      ag.mon.mbx_mon_to_sb  = mbx_mon_to_sb;
      ag.mon.mbx_mon_to_cov = mbx_mon_to_cov;

    endfunction : build

    // =========================
    // RUN PHASE
    // =========================
    task run();
      fork
        ag.run();   // Generator + driver + monitor
        sb.run();   // Scoreboard
        cov.run();  // Coverage
      join_none

      #2100;
      // Final report of error
      sb.report();
      
    endtask : run

  endclass : ALU_env

endpackage : ALU_env_pkg
