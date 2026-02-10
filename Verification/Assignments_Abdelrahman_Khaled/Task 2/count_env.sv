`timescale 1ns / 1ps

package count_env_pkg;

  // Import required packages
  import count_transaction_pkg::*;
  import count_agent_pkg::*;
  import count_coverage_pkg::*;
  import count_scoreboard_pkg::*;

  // ---------------------------------------------------------------------------
  // count_env class
  // Top-level environment that connects agent, scoreboard, and coverage
  // ---------------------------------------------------------------------------
  class count_env;

    // =========================
    // COMPONENTS
    // =========================
    count_agent       ag;   // Agent (gen + drv + mon)
    count_scoreboard  sb;   // Scoreboard
    count_coverage    cov;  // Coverage

    // =========================
    // MAILBOXES
    // =========================
    mailbox #(count_transaction) mbx_mon_to_sb;
    mailbox #(count_transaction) mbx_mon_to_cov;

    // =========================
    // VIRTUAL INTERFACES
    // =========================
    virtual count_if.tb      vif_tb;   // Driver interface
    virtual count_if.monitor vif_mon;  // Monitor interface

    // =========================
    // CONSTRUCTOR (empty)
    // =========================
    function new();
      // Create environment-level mailboxes
      mbx_mon_to_sb  = new();
      mbx_mon_to_cov = new();
    endfunction : new

    // =========================
    // CONNECT FUNCTION
    // =========================
    function void connect(
      virtual count_if.tb      vif_tb,
      virtual count_if.monitor vif_mon
    );
      // Save virtual interfaces
      this.vif_tb  = vif_tb;
      this.vif_mon = vif_mon;

      // -------------------------
      // Create & connect agent
      // -------------------------
      ag = new();
      ag.connect(vif_tb, vif_mon);

      // -------------------------
      // Create scoreboard & coverage
      // -------------------------
      sb  = new();
      cov = new();

      // Connect mailboxes
      sb.connect(mbx_mon_to_sb);
      cov.connect(mbx_mon_to_cov);

      // Connect agent monitor outputs
      ag.mon.mbx_mon_to_sb  = mbx_mon_to_sb;
      ag.mon.mbx_mon_to_cov = mbx_mon_to_cov;

    endfunction : connect

    // =========================
    // RUN PHASE
    // =========================
    task run();
      fork
        ag.run();   // Generator + Driver + Monitor
        sb.run();   // Scoreboard
        cov.run();  // Coverage
      join_none

      // Let simulation run
      #24000;

      // Final report
      sb.report();
    endtask : run

  endclass : count_env

endpackage : count_env_pkg
