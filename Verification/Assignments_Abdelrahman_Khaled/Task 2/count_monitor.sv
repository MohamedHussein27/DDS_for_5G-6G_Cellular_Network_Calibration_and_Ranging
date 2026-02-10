// -----------------------------------------------------------------------------
// Package: count_monitor_pkg
// Description:
//   This package defines a monitor for the count FSM. The monitor observes DUT
//   signals through a virtual interface, packs them into a transaction object,
//   and sends the transaction to the scoreboard and coverage components.
// -----------------------------------------------------------------------------
`timescale 1ns / 1ps

package count_monitor_pkg;

  import count_transaction_pkg::*;

  // ---------------------------------------------------------------------------
  // Count Monitor class
  // Observes DUT signals and creates transactions for verification
  // ---------------------------------------------------------------------------
  class count_monitor;

    // -----------------------------
    // Transaction object to hold observed values
    // -----------------------------
    count_transaction tr_mon; 

    // -----------------------------
    // Virtual interface handle (monitor modport)
    // -----------------------------
    virtual count_if.monitor vif;

    // -----------------------------
    // Mailboxes to communicate with scoreboard and coverage
    // -----------------------------
    mailbox #(count_transaction) mbx_mon_to_sb;
    mailbox #(count_transaction) mbx_mon_to_cov;

    // -----------------------------
    // Constructor: initialize transaction
    // -----------------------------
    function new();
      tr_mon = new();
    endfunction : new

    // -----------------------------
    // Connect function: attach interface and mailboxes
    // -----------------------------
    function void connect(
      virtual count_if.monitor vif,
      mailbox #(count_transaction) mbx_mon_to_sb,
      mailbox #(count_transaction) mbx_mon_to_cov
    );
      this.vif            = vif;
      this.mbx_mon_to_sb  = mbx_mon_to_sb;
      this.mbx_mon_to_cov = mbx_mon_to_cov;
    endfunction : connect

    // -----------------------------
    // Main monitoring task
    // -----------------------------
    task run();
      forever begin
        // Wait for positive edge of DUT clock
        @(posedge vif.clk);
        #0.01; // small delay to sample stable signals

        // -----------------------------
        // Sample DUT outputs into transaction object
        // -----------------------------
        tr_mon.rst_n       = vif.rst_n;
        tr_mon.start       = vif.start;
        tr_mon.flag        = vif.flag;
        tr_mon.wait_timer  = vif.wait_timer;
        tr_mon.busy        = vif.busy;
        tr_mon.count_value = vif.count_value;

        // -----------------------------
        // Send transaction copies to consumers
        // -----------------------------
        mbx_mon_to_sb.put(tr_mon);
        mbx_mon_to_cov.put(tr_mon);
      end
    endtask : run

  endclass : count_monitor

endpackage : count_monitor_pkg
