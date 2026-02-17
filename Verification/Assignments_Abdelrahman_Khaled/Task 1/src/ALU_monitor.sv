package ALU_monitor_pkg;

  // Import required packages
  import ALU_transaction_pkg::*;

  // ALU Monitor class
  // Passively observes DUT signals and forwards transactions
  class ALU_monitor;

    // Virtual interface handle (monitor modport)
    virtual ALU_if.monitor vif;

    // Mailboxes to scoreboard and coverage
    mailbox #(ALU_transaction) mbx_mon_to_sb;
    mailbox #(ALU_transaction) mbx_mon_to_cov;

    // Constructor
    function new(
      virtual ALU_if.monitor vif,
      mailbox #(ALU_transaction) mbx_mon_to_sb,
      mailbox #(ALU_transaction) mbx_mon_to_cov
    );
      this.vif            = vif;
      this.mbx_mon_to_sb  = mbx_mon_to_sb;
      this.mbx_mon_to_cov = mbx_mon_to_cov;
    endfunction : new

    // Main monitoring task
    task run();
      ALU_transaction tr_mon; // local transaction (new each cycle)

      forever begin
          
        // Synchronize sampling
        @(posedge vif.clk);
        #1;

        // Create a fresh transaction for this sample
        tr_mon = new();

        // Sample ALU inputs and outputs
        tr_mon.a   = vif.a;
        tr_mon.b   = vif.b;
        tr_mon.op  = vif.op;
        tr_mon.out = vif.out;
        tr_mon.c   = vif.c;

        // Send independent copies to consumers
        mbx_mon_to_sb.put(tr_mon);
        mbx_mon_to_cov.put(tr_mon);

      end
    endtask : run

  endclass : ALU_monitor

endpackage : ALU_monitor_pkg
