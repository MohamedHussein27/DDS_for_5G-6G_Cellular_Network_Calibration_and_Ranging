
package counter_sequence_pkg;
  import counter_seq_item_pkg::*;

  class counter_sequence;

    mailbox #(counter_seq_item) seq2drv_mb;
    counter_seq_item item;

    function new();
    endfunction

    function void connect(mailbox #(counter_seq_item) mb);
      seq2drv_mb = mb;
    endfunction

    // ---------------------------------------------------
    // SEND TASK (simple, clean, your style)
    // ---------------------------------------------------
    task send(
      logic        rst_n,
      logic        start,
      logic        flag,
      logic [15:0] wait_timer
    );
      item = new();
      item.rst_n      = rst_n;
      item.start      = start;
      item.flag       = flag;
      item.wait_timer = wait_timer;
      seq2drv_mb.put(item);
    endtask

    // ---------------------------------------------------
    // MAIN SEQUENCE
    // ---------------------------------------------------
    task run();

      
      // =====================================================
      // SEQ1 : Count normally until MAX (31)
      // =====================================================
      $display("=== SEQ1: Count to max ===");

      // start ↑ (ONE cycle only)
      send(1, 1, 0, 3);

      // drop start
      send(1, 0, 0, 3);

      // allow counter to reach max
      repeat (150)
        send(1, 0, 0, 3);

     

     // =====================================================
      // SEQ2 : Assert FLAG exactly in the CHECK state
      // =====================================================
      $display("=== SEQ2: Flag in CHECK state ===");
      
      // 1. Send START (Cycle 0: Scoreboard goes IDLE -> WAIT)
      send(1, 1, 0, 3); 

      // 2. Wait exactly 'wait_timer' cycles (Cycles 1, 2, 3)
      // Scoreboard timer counts down: 3, then 2, then 1.
      repeat (3) begin
        send(1, 0, 0, 3);
      end

      // 3. Assert FLAG (Cycle 4)
      // Scoreboard timer finished last cycle, so it enters SB_CHECK right NOW.
      // It will see the flag=1 and go straight back to IDLE.
      send(1, 0, 1, 3);

      // 4. Clear flag and Idle out
      send(1, 0, 0, 3);
      repeat (5) begin
        send(1, 0, 0, 3);
      end


      // -----------------------------------------------------
      // TEST 3: Zero Wait Timer (Boundary Condition)
      // Spec says "programmable number"
      // -----------------------------------------------------
      $display("--- TEST 3: Zero Wait Timer ---");
      send(1, 1, 0, 0); // Pulse start with timer=0
      send(1, 0, 0, 0);
      repeat (40) send(1, 0, 0, 0); // Should count very fast


      // -----------------------------------------------------
      // TEST 4:
      // Asserts flag at the exact cycle the timer expires.
      // -----------------------------------------------------
      $display("--- TEST 4: Flag precisely on CHECK state ---");
      send(1, 1, 0, 3); // Start (Cycle 0)
      repeat (3) send(1, 0, 0, 3); // Wait countdown 
      send(1, 0, 1, 3); // FLAG on Cycle 4 (CHECK state)
      send(1, 0, 0, 3);
      repeat (5) send(1, 0, 0, 3);


      
      // TEST 5: Start Glitch (Re-asserting while busy)
      // Does the FSM ignore a new start, or does it break
      
      $display("--- TEST 5: Start pulse while already busy ---");
      send(1, 1, 0, 5); // Normal start
      send(1, 0, 0, 5);
      repeat (4) send(1, 0, 0, 5); 
      send(1, 1, 0, 5); // ROGUE START PULSE!
      send(1, 0, 0, 5);
      repeat (10) send(1, 0, 0, 5);
      send(1, 0, 1, 5); // Cleanup via flag
      send(1, 0, 0, 5);


      // -----------------------------------------------------
      // TEST 6: Start & Flag Simultaneously (Priority Check)
      // Tests which signal wins if both arrive on the same edge.
      // -----------------------------------------------------
      $display("--- TEST 6: Simultaneous Start and Flag ---");
      send(1, 1, 1, 3); // Start=1 AND Flag=1
      send(1, 0, 0, 3);
      repeat (5) send(1, 0, 0, 3);


      
      // TEST 7: Asynchronous Reset Mid-Flight
      // Drops reset to 0 while the FSM is actively counting.
      
      $display("--- TEST 7: Asynchronous Reset Recovery ---");
      send(1, 1, 0, 3);
      send(1, 0, 0, 3);
      repeat (8) send(1, 0, 0, 3);
      send(0, 0, 0, 3); // DROP RESET!
      send(1, 0, 0, 3); // Bring reset back high
      repeat (5) send(1, 0, 0, 3);


      // -----------------------------------------------------
      // TEST 8: Changing wait_timer while busy
      // Checks if the hardware properly latches the timer input.
      // -----------------------------------------------------
      $display("--- TEST 8: Mutating wait_timer mid-count ---");
      send(1, 1, 0, 5); // Start with timer = 5
      send(1, 0, 0, 5);
      repeat (6) send(1, 0, 0, 5);
      send(1, 0, 0, 1); // Mutate timer to 1!
      repeat (10) send(1, 0, 0, 1);
      send(1, 0, 1, 1); // Cleanup  flag
      send(1, 0, 0, 1);


      // -----------------------------------------------------
      // TEST 9: Back-to-Back Starts
      // Asserts start immediately after a sequence ends.
      // -----------------------------------------------------
      $display("--- TEST 9: Back-to-Back Sequences ---");
      send(1, 1, 0, 2); 
      send(1, 0, 0, 2);
      repeat (5) send(1, 0, 0, 2);
      send(1, 0, 1, 2); // Stop first sequence
      send(1, 1, 0, 3); // IMMEDIATELY start a new sequence
      send(1, 0, 0, 3);
      repeat (10) send(1, 0, 0, 3);
      send(1, 0, 1, 3); // Cleanup
      send(1, 0, 0, 3);


      // -----------------------------------------------------
      // TEST 10: Larger wait_timer (Datapath Check)
      // Ensures the 16-bit timer doesn't fail on double-digit numbers.
      // -----------------------------------------------------
      $display("--- TEST 10: Larger wait_timer (10) ---");
      send(1, 1, 0, 10);
      send(1, 0, 0, 10);
      repeat (25) send(1, 0, 0, 10); // Let it count twice
      send(1, 0, 1, 10); // Cleanup
      send(1, 0, 0, 10);
      repeat (5) send(1, 0, 0, 10);

     
      $display(" ALL DIRECTED TESTS COMPLETE ");
      

     
    endtask

  endclass
endpackage
