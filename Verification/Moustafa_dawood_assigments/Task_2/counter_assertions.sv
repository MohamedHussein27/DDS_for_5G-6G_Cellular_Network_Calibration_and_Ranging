module counter_assertions (counter_bfm bfm);

  // ---------------------------------------------------
  // Default clocking & Global Disable
  // ---------------------------------------------------
  default clocking cb @(posedge bfm.clk); endclocking
  
  // Disable all assertions if reset is active (0) to prevent false failures
  default disable iff (!bfm.rst_n);

  // ---------------------------------------------------
  // 1. ASYNC RESET (Check Combinational Output)
  // When rst_n goes low, busy and count_value must immediately be 0.
  // We don't use the clock (->) because it's asynchronous.
  // ---------------------------------------------------
  property p_async_reset;
    @(bfm.rst_n) (!bfm.rst_n) |-> (bfm.busy == 0 && bfm.count_value == 0);
  endproperty

  assert property (p_async_reset)
    else $error("ASSERTION FAIL: Reset did not immediately clear busy and count");
  cover property (p_async_reset);


  // ---------------------------------------------------
  // 2. START TO BUSY (1-Cycle Delay)
  // A rising edge on start (when idle) MUST trigger busy on the NEXT cycle.
  // ---------------------------------------------------
  property p_start_causes_busy;
    ($rose(bfm.start) && !bfm.busy) |=> bfm.busy;
  endproperty

  assert property (p_start_causes_busy)
    else $error("ASSERTION FAIL: Rising edge on start did not assert busy");
  cover property (p_start_causes_busy);


  // ---------------------------------------------------
  // 3. FLAG ABORT (Immediate Stop)
  // If flag is asserted while busy, the very next cycle busy MUST drop.
  // ---------------------------------------------------
  property p_flag_causes_immediate_stop;
    (bfm.busy && bfm.flag) |=> !bfm.busy;
  endproperty

  assert property (p_flag_causes_immediate_stop)
    else $error("ASSERTION FAIL: Flag did not immediately abort busy");
  cover property (p_flag_causes_immediate_stop);


  // ---------------------------------------------------
  // 4. MAX COUNT AUTO-STOP
  // If the counter hits 31, it MUST drop busy on the next clock cycle.
  // ---------------------------------------------------
  property p_stop_at_max_count;
    (bfm.count_value == 5'd31 && bfm.busy) |=> !bfm.busy;
  endproperty

  assert property (p_stop_at_max_count)
    else $error("ASSERTION FAIL: FSM did not drop busy after hitting max count (31)");
  cover property (p_stop_at_max_count);


  // ---------------------------------------------------
  // 5. HOLD VALUE WHILE IDLE
  // If the FSM is not busy, the count_value must NOT change, 
  // UNLESS a new start edge occurs.
  // ---------------------------------------------------
  property p_hold_value_when_idle;
    (!bfm.busy && !$rose(bfm.start)) |=> (bfm.count_value == $past(bfm.count_value));
  endproperty

  assert property (p_hold_value_when_idle)
    else $error("ASSERTION FAIL: count_value changed while FSM was idle");
  cover property (p_hold_value_when_idle);


  // ---------------------------------------------------
  // 6. CLEAR COUNT ON NEW START
  // A valid start edge must clear the counter to 0 on the next cycle.
  // ---------------------------------------------------
  property p_clear_count_on_start;
    ($rose(bfm.start) && !bfm.busy) |=> (bfm.count_value == 0);
  endproperty

  assert property (p_clear_count_on_start)
    else $error("ASSERTION FAIL: New start edge did not clear count_value to 0");
  cover property (p_clear_count_on_start);


  // ---------------------------------------------------
  // 7. NO OVERFLOW
  // The 5-bit counter should mathematically never exceed 31, 
  // but it's a good safety check against weird RTL rollover bugs.
  // ---------------------------------------------------
  property p_no_overflow;
    bfm.count_value <= 5'd31;
  endproperty

  assert property (p_no_overflow)
    else $error("ASSERTION FAIL: count_value exceeded 31");
  cover property (p_no_overflow);

endmodule