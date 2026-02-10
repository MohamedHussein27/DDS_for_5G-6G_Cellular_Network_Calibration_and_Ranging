/*
package counter_scoreboard_pkg;

  import counter_seq_item_pkg::*;
  import counter_globals_pkg::*;

  class counter_scoreboard;

    mailbox #(counter_seq_item) mon2sb_mb;

    // ----------------------------------------
    // Constructor
    // ----------------------------------------
    function new();
      sb_state     = SB_IDLE;
      sb_count     = 0;
      sb_timer     = 0;
      sb_busy      = 0;
      sb_start     = 0;
      sb_flag      = 0;
      sb_wait_timer= 0;
      prev_start   = 0;
      error_count  = 0;
    endfunction

    // ----------------------------------------
    // Main scoreboard loop
    // ----------------------------------------
    task run();
      counter_seq_item item;
      logic start_rise;

      forever begin
        // Get monitored transaction
        mon2sb_mb.get(item);

        // ------------------------------------
        // SAMPLE DUT inputs (ONLY place item.* is allowed)
        // ------------------------------------
        sb_start       = item.start;
        sb_flag        = item.flag;
        sb_wait_timer  = item.wait_timer;

        // Detect rising edge of start
        start_rise = sb_start && !prev_start;
        prev_start = sb_start;

        // ------------------------------------
        // Reference FSM (sb_* ONLY)
        // ------------------------------------
        if (!item.rst_n) begin
          sb_state = SB_IDLE;
          sb_count = 0;
          sb_timer = 0;
        end
        else begin
          case (sb_state)

            // ---------------- IDLE ----------------
            SB_IDLE: begin
              sb_count = 0;
              if (start_rise) begin
                sb_timer = sb_wait_timer;
                sb_state = SB_WAIT;
              end
            end

            // ---------------- WAIT ----------------
            SB_WAIT: begin
              if (sb_timer > 0)
                sb_timer--;
              else
                sb_state = SB_CHECK;
            end

            // ---------------- CHECK ----------------
            SB_CHECK: begin
              if (sb_flag || sb_count == 5'd31) begin
                sb_state = SB_STOP;
              end
              else begin
                sb_count++;
                sb_timer = sb_wait_timer;
                sb_state = SB_WAIT;
              end
            end

            // ---------------- STOP ----------------
            SB_STOP: begin
              sb_state = SB_IDLE; // one cycle only
            end

          endcase
        end

        // ------------------------------------
        // Busy definition
        // ------------------------------------
        sb_busy = (sb_state != SB_IDLE);

        // ------------------------------------
        // DUT vs Reference comparison
        // ------------------------------------
        if (item.busy !== sb_busy) begin
          error_count++;
          $display("ERROR BUSY  EXP=%0b GOT=%0b  errors=%0d",
                   sb_busy, item.busy, error_count);
        end

        if (item.count_value !== sb_count) begin
          error_count++;
          $display("ERROR COUNT EXP=%0d GOT=%0d errors=%0d",
                   sb_count, item.count_value, error_count);
        end

        if ((item.busy === sb_busy) &&
            (item.count_value === sb_count)) begin
          $display("PASS busy=%0b count=%0d",
                   item.busy, item.count_value);
        end
      end
    endtask

  endclass

endpackage
*/
/*
package counter_scoreboard_pkg;

  import counter_seq_item_pkg::*;
  import counter_globals_pkg::*;

  class counter_scoreboard;

    mailbox #(counter_seq_item) mon2sb_mb;

    // ----------------------------------------
    // Constructor
    // ----------------------------------------
    function new();
      sb_state     = SB_IDLE;
      sb_count     = 0;
      sb_timer     = 0;
      sb_busy      = 0;
      sb_start     = 0;
      sb_flag      = 0;
      sb_wait_timer= 0;
      prev_start   = 0;
      error_count  = 0;
    endfunction

    // ----------------------------------------
    // Main scoreboard loop
    // ----------------------------------------
    task run();
      counter_seq_item item;
      logic start_rise;

      forever begin
        mon2sb_mb.get(item);

        // ------------------------------------
        // SAMPLE INPUTS
        // ------------------------------------
        sb_start       = item.start;
        sb_flag        = item.flag;
        sb_wait_timer  = item.wait_timer;

        // Rising Edge Detection
        start_rise = sb_start && !prev_start;
        prev_start = sb_start;

        // ------------------------------------
        // Reference FSM
        // ------------------------------------
        if (!item.rst_n) begin
          sb_state = SB_IDLE;
          sb_count = 0;
          sb_timer = 0;
        end
        else begin
          case (sb_state)

            // ---------------- IDLE ----------------
            SB_IDLE: begin
              sb_count = 0;
              if (start_rise) begin
                sb_timer = sb_wait_timer;
                sb_state = SB_WAIT;
              end
            end

            // ---------------- WAIT ----------------
            SB_WAIT: begin
              
              // FIX: Allow timer to reach 0 before acting.
              // This consumes exactly 4 clock cycles (for wait_timer=3),
              // aligning perfectly with the DUT.
              if (sb_timer > 0) begin
                sb_timer--;
              end
              else begin
                // Timer is 0. Perform check and increment.
                
                if (sb_flag || sb_count == 5'd31) begin
                   sb_state = SB_IDLE; 
                end
                else begin
                   sb_count++;
                   sb_timer = sb_wait_timer; 
                end
              end
            end

          endcase
        end

        // ------------------------------------
        // Busy Logic
        // ------------------------------------
        sb_busy = (sb_state != SB_IDLE);

        // ------------------------------------
        // Comparison Logic
        // ------------------------------------
        if (item.busy !== sb_busy) begin
          error_count++;
          $display("ERROR BUSY  EXP=%0b GOT=%0b  errors=%0d", 
                   sb_busy, item.busy, error_count);
        end

        if (item.count_value !== sb_count) begin
          error_count++;
          $display("ERROR COUNT EXP=%0d GOT=%0d errors=%0d", 
                   sb_count, item.count_value, error_count);
        end

        if ((item.busy === sb_busy) && (item.count_value === sb_count)) begin
          $display("PASS busy=%0b count=%0d", item.busy, item.count_value);
        end

      end
    endtask

  endclass

endpackage
*/
package counter_scoreboard_pkg;

  import counter_seq_item_pkg::*;
  import counter_globals_pkg::*;

  class counter_scoreboard;

    mailbox #(counter_seq_item) mon2sb_mb;

    // ----------------------------------------
    // FSM States (Clean & Canonical)
    // ----------------------------------------
    
    

    
    function new();
      sb_state     = SB_IDLE;
      sb_count     = 0;
      sb_timer     = 0;
      sb_busy      = 0;
      sb_start     = 0;
      sb_flag      = 0;
      sb_wait_timer= 0;
      prev_start   = 0;
      error_count  = 0;
    endfunction

    task run();
      counter_seq_item item;
      logic start_rise;

      forever begin
        mon2sb_mb.get(item);


        if (item.busy !== sb_busy) begin
          error_count++;
          $display("ERROR BUSY  EXP=%0b GOT=%0b  errors=%0d time=%0t", 
                   sb_busy, item.busy, error_count, $time);
        end

        if (item.count_value !== sb_count) begin
          error_count++;
          $display("ERROR COUNT EXP=%0d GOT=%0d errors=%0d time=%0t", 
                   sb_count, item.count_value, error_count, $time);
        end
        
        if ((item.busy === sb_busy) && (item.count_value === sb_count)) begin
             $display("PASS busy=%0b count=%0d", item.busy, item.count_value);
        end
        

        // Sample Inputs
        sb_start       = item.start;
        sb_flag        = item.flag;
        sb_wait_timer  = item.wait_timer;
        

        // Rising Edge Detection
        start_rise = sb_start && !prev_start;
        prev_start = sb_start;

        if (!item.rst_n) begin
          sb_state = SB_IDLE;
          sb_count = 0;
          sb_timer = 0;
        end
        else begin
          case (sb_state)

            // ---------------- IDLE ----------------
            // ---------------- IDLE ----------------
// ---------------- IDLE ----------------
            SB_IDLE: begin
              sb_busy = 0; // Busy is low while idling

              if (start_rise) begin
                // 1. CLEAR THE COUNT HERE
                // The count only goes to 0 when a brand new sequence starts
                sb_count = 0;
                
                // 2. Load the full timer first
                sb_timer = sb_wait_timer; 
                
                // 3. Go to WAIT
                sb_state = SB_WAIT;
              end
            end

            // ---------------- WAIT ----------------
            SB_WAIT: begin
              
              sb_busy = 1; // Busy is high in WAIT and CHECK states
              if (sb_timer > 1) begin
                sb_timer--;
              end
              else begin
                // Timer finished, go to Check
                sb_state = SB_CHECK;
              end
            end

            // ---------------- CHECK ----------------
            // Takes 1 Cycle. Handles Logic & Reloading.
           // ---------------- CHECK ----------------
            // Takes 1 Cycle. Handles Logic & Reloading.
            SB_CHECK: begin
              
              // 1. Check Flag first
              if (sb_flag) begin
                sb_state = SB_IDLE;
              end
              else begin
                
                
                    sb_count++;
                
                // 3. IMMEDIATELY check if we hit the max (31)
                if (sb_count == 5'd31) begin
                  sb_state = SB_IDLE; // Drop to IDLE now. Busy goes low next cycle.
                end
                else begin
                  // 4. Only reload and wait if we haven't reached 31
                  sb_timer = sb_wait_timer - 1;
                  sb_state = SB_WAIT;
                end
                
              end
            end

          endcase
        end

        // ------------------------------------
        // Busy Logic
        // ------------------------------------
        //sb_busy = (sb_state != SB_IDLE);

        // ------------------------------------
        // Comparison
        // ------------------------------------
      

      end
    endtask

  endclass

endpackage
