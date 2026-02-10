package count_scoreboard_pkg;

  import count_transaction_pkg::*;

  typedef class IdleState;
  typedef class ReadyState;
  typedef class CountingState;

  ///////////////////////////////////////////////////////////////////
  // Reference signals
  logic       busy_ref        = 0;
  logic [4:0] count_value_ref = 0;
  int         timer_counter   = 0;
  int         flag_counter    = 0;

  // FSM states enumeration
  typedef enum logic [1:0] {
    IDLE,
    READY,
    COUNTING
  } STATES;

  STATES current_state;
  STATES next_state;

  ///////////////////////////////////////////////////////////////////
  // Base State class
  class State;

    virtual function State transition(
      logic rst_n,
      logic start,
      logic flag,
      logic [15:0] wait_timer
    );
      return this;
    endfunction

    function State reset_triggered();
      IdleState next_state;
      next_state = new();
      return next_state;
    endfunction

    function State go_to_idle_state();
      IdleState next_state = new();
      return next_state;
    endfunction

    function State go_to_ready_state();
      ReadyState next_state = new();
      return next_state;
    endfunction

    function State go_to_count_state();
      CountingState next_state = new();
      return next_state;
    endfunction

  endclass : State

  ///////////////////////////////////////////////////////////////////
  // IDLE STATE
  ///////////////////////////////////////////////////////////////////
  class IdleState extends State;

    function State transition(
      logic rst_n,
      logic start,
      logic flag,
      logic [15:0] wait_timer
    );
      current_state   = IDLE;
      busy_ref        = 0;
      timer_counter   = 0;
      flag_counter    = 0;

      if (!rst_n) begin
        count_value_ref = 0;
        next_state = IDLE;
        return reset_triggered();
      end

      if (start) begin
        next_state = READY;
        return go_to_ready_state();   // IDLE → READY
      end

      return this;
    endfunction

  endclass : IdleState

  ///////////////////////////////////////////////////////////////////
  // READY STATE
  ///////////////////////////////////////////////////////////////////
  class ReadyState extends State;

    function State transition(
      logic rst_n,
      logic start,
      logic flag,
      logic [15:0] wait_timer
    );
      current_state   = READY;
      busy_ref        = 1;   // busy asserted
      count_value_ref = 0;   // count cleared
      timer_counter   = 0;
      flag_counter    = 0;

      if (!rst_n) begin
        busy_ref        = 0;
        next_state = IDLE;
        return reset_triggered();
        end

      next_state = COUNTING;
      // READY lasts exactly one cycle
      return go_to_count_state();   // READY → COUNTING
    endfunction

  endclass : ReadyState

  ///////////////////////////////////////////////////////////////////
  // COUNTING STATE
  ///////////////////////////////////////////////////////////////////
  class CountingState extends State;

    function State transition(
      logic rst_n,
      logic start,
      logic flag,
      logic [15:0] wait_timer
    );

      current_state = COUNTING;
      busy_ref      = 1;

      // ------------------------------------------------------------
      // Counters update
      // ------------------------------------------------------------
      timer_counter++;
      flag_counter++;

      // ------------------------------------------------------------
      // Count increment logic
      // ------------------------------------------------------------
      if (timer_counter == wait_timer) begin
        timer_counter = 0;

        if (count_value_ref < 5'd31)
          count_value_ref++;
      end

      // Align flag counter early
      if (timer_counter == 1)
        flag_counter = timer_counter;

      // ------------------------------------------------------------
      // STOP CONDITIONS
      // ------------------------------------------------------------

      // Reset or saturation
      if (!rst_n) begin
        timer_counter = 0;
        flag_counter  = 0;
        busy_ref        = 0;
        count_value_ref = 0;
        next_state = IDLE;
        return reset_triggered();
      end

      if (count_value_ref == 5'd31) begin
        timer_counter = 0;
        flag_counter  = 0;
        next_state = IDLE;
        return go_to_idle_state();
      end

      // Flag at interval boundary
      if ((flag_counter == wait_timer) && flag) begin
        timer_counter = 0;
        flag_counter  = 0;
        count_value_ref--;  // compensate
        next_state = IDLE;
        return go_to_idle_state();
      end

      return this;

    endfunction

  endclass : CountingState

  ///////////////////////////////////////////////////////////////////
  // Scoreboard class
  class count_scoreboard;

    mailbox #(count_transaction) mbx_mon_to_sb;
    count_transaction tr_sb;

    int Error_count   = 0;
    int Correct_count = 0;

    IdleState starting_state;
    State current_state_obj;
    State next_state_obj;

    // Constructor function
    function new();
      starting_state = new();
      this.current_state_obj = starting_state;
      this.next_state_obj    = current_state_obj;
    endfunction
    
    // -----------------------------
    // Connect function: attach mailbox and initialize state
    // -----------------------------
    function void connect(mailbox #(count_transaction) mbx);
      this.mbx_mon_to_sb  = mbx;
    endfunction

    // Compare DUT transaction to reference model
    function void count_check(count_transaction tr_DUT);
      // Update reference model
      next_state_obj = current_state_obj.transition(
        tr_DUT.rst_n, tr_DUT.start, tr_DUT.flag, tr_DUT.wait_timer
      );
      current_state_obj = next_state_obj;

      // Compare DUT vs reference
      if (tr_DUT.busy !== busy_ref || tr_DUT.count_value !== count_value_ref) begin
        Error_count++;
        $display("ERROR at time %0t: DUT={busy=%b count=%0d} REF={busy=%b count=%0d}",
                  $time, tr_DUT.busy, tr_DUT.count_value, busy_ref, count_value_ref);
      end else begin
        Correct_count++;
      end
    endfunction

    // Main run task
    task run();
      forever begin
        mbx_mon_to_sb.get(tr_sb);
        count_check(tr_sb);
      end
    endtask

    // Final report
    task report();
      $display("Simulation summary: Correct=%0d Errors=%0d", Correct_count, Error_count);
    endtask

  endclass : count_scoreboard

endpackage : count_scoreboard_pkg
