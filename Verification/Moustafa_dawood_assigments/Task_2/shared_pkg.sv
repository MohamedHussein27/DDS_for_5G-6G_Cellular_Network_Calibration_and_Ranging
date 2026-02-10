package counter_globals_pkg;

  // ----------------------------------------
  // Reference FSM states
  // ----------------------------------------
  typedef enum logic [1:0] {
    SB_IDLE,
    SB_WAIT,
    SB_CHECK,
    SB_STOP
  } sb_state_t;

  // ----------------------------------------
  // Scoreboard-visible sampled inputs
  // ----------------------------------------
  logic       sb_start;
  logic       sb_flag;
  logic [4:0] sb_wait_timer;

  // ----------------------------------------
  // Reference model internal state
  // ----------------------------------------
  sb_state_t  sb_state;
  logic [4:0] sb_count;
  logic [4:0] sb_timer;
  logic       sb_busy;
  logic       prev_start;

  // ----------------------------------------
  // Error counter
  // ----------------------------------------
  int error_count;

endpackage
