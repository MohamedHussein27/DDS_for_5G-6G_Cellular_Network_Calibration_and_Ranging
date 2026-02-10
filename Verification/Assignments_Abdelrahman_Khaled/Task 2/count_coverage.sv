package count_coverage_pkg;

  import count_transaction_pkg::*;

  // ---------------------------------------------------------------------------
  // count_coverage class
  // Implements functional coverage for the FSM counter transactions
  // ---------------------------------------------------------------------------
  class count_coverage;

    // Transaction handle for covergroup sampling
    count_transaction tr_cvg;

    // Mailbox connected to monitor
    mailbox #(count_transaction) mbx_mon_to_cov;

    // Covergroup declaration
    covergroup count_cg;

      // Coverpoint for start signal
      coverpoint tr_cvg.start {
        bins start_0 = {0};
        bins start_1 = {1};
      }

      // Coverpoint for flag signal
      coverpoint tr_cvg.flag {
        bins flag_0 = {0};
        bins flag_1 = {1};
      }

      // Coverpoint for wait_timer (example bins)
      coverpoint tr_cvg.wait_timer {
      bins all_values = {[0:65535]};
}


      // Coverpoint for busy
      coverpoint tr_cvg.busy {
        bins busy_0 = {0};
        bins busy_1 = {1};
      }

      // Coverpoint for count_value
      coverpoint tr_cvg.count_value {
        bins cnt_bins[] = {[0:31]};
      }

    endgroup : count_cg

    // Constructor: initialize transaction and covergroup
    function new();
      tr_cvg = new();
      count_cg = new(); // bind to tr_cvg automatically
    endfunction : new

    // -----------------------------
    // Connect function: assign mailbox
    // -----------------------------
    function void connect(mailbox #(count_transaction) mbx_mon_to_cov);
      this.mbx_mon_to_cov = mbx_mon_to_cov;
    endfunction : connect

    // -----------------------------
    // Sample a transaction into the covergroup
    // -----------------------------
    function void sample_data(count_transaction tr_DUT);
      tr_cvg.start       = tr_DUT.start;
      tr_cvg.flag        = tr_DUT.flag;
      tr_cvg.wait_timer  = tr_DUT.wait_timer;
      tr_cvg.busy        = tr_DUT.busy;
      tr_cvg.count_value = tr_DUT.count_value;

      count_cg.sample();
    endfunction : sample_data

    // -----------------------------
    // Main run task: continuously sample transactions from monitor
    // -----------------------------
    task run();
      forever begin
        
        mbx_mon_to_cov.get(tr_cvg); // Get new transaction from monitor
        sample_data(tr_cvg);        // Sample it into covergroup
      end
    endtask : run

  endclass : count_coverage

endpackage : count_coverage_pkg
