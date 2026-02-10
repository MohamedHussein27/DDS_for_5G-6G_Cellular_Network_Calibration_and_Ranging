// -----------------------------------------------------------------------------
// Package: count_transaction_pkg
// Description:
//   This package defines the transaction class for the count FSM. Each instance
//   of count_transaction represents one stimulus item applied to the count FSM.
//   It contains randomized input signals and observed DUT outputs for verification.
// -----------------------------------------------------------------------------
package count_transaction_pkg;

  // ------------------------------------------------------------
  // Count Transaction
  // Represents one stimulus item applied to the count FSM
  // ------------------------------------------------------------
  class count_transaction;

    // -----------------------------
    // Control / configuration inputs
    // -----------------------------

    // Active-low asynchronous reset
    rand logic rst_n;

    // Start pulse: triggers counting, typically asserted for one cycle
    rand logic start;

    // Flag that may terminate counting early (used by FSM logic)
    rand logic flag;

    // Programmable wait interval (number of clock cycles to wait between counts)
    rand logic [15:0] wait_timer;

    // -----------------------------
    // Constraints for input randomization
    // -----------------------------

    // rst_n: mostly 1 (no reset), rarely 0 (reset)
    constraint reset {
      rst_n dist { 0 := 2, 1 := 98 }; // 2% chance reset active, 98% chance inactive
    }

    // start: mostly asserted, sometimes not asserted
    constraint start_count {
      start dist { 0 := 40, 1 := 60 }; // 60% chance start pulse, 40% no pulse
    }

    // flag: mostly 0 , rarely 1 
    constraint flag_count {
      flag dist { 0 := 85, 1 := 15 }; // 15% chance reset active, 85% chance inactive
    }

    // wait_timer: mostly 1000 cycles, occasionally slightly lower or higher
    constraint register {
        wait_timer == 5; // always 1000
    }


    // -----------------------------
    // DUT outputs (observed, not randomized)
    // -----------------------------

    // Indicates that FSM is actively counting
    logic busy;

    // Count value produced by the FSM
    logic [4:0] count_value;

  endclass : count_transaction

endpackage : count_transaction_pkg
