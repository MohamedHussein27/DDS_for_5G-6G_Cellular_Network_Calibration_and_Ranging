package count_generator_pkg;
  import count_transaction_pkg::*;

  class count_generator;
    
    count_transaction tr_gen;
    mailbox #(count_transaction) gen_to_dr;

    function new();
    endfunction

    function void connect(mailbox #(count_transaction) gen_to_dr);
      this.gen_to_dr = gen_to_dr;
    endfunction

    task send(
      logic rst_n,
      logic start,
      logic flag,
      logic [15:0] wait_timer
    );
      tr_gen = new();
      tr_gen.rst_n      = rst_n;
      tr_gen.start      = start;
      tr_gen.flag       = flag;
      tr_gen.wait_timer = wait_timer;
      gen_to_dr.put(tr_gen);
    endtask

    task run();

      // =====================================================
      // 1. RESET PHASE
      // =====================================================
      repeat (5)
        send(0, 0, 0, 0);

      // Release reset
      send(1, 0, 0, 0);

      // =====================================================
      // 2. DIRECTED FSM TRANSITIONS
      // =====================================================

      // Start counting
      send(1, 1, 0, 5);

      // Let counter run
      repeat (10)
        send(1, 0, 0, 5);

      // Stop via flag
      send(1, 0, 1, 5);

      // Idle cycles
      repeat (5)
        send(1, 0, 0, 5);

      // =====================================================
      // 3. COUNTER MAX TEST
      // =====================================================
      send(1, 1, 0, 1);      // fast count
      repeat (40)
        send(1, 0, 0, 1);    // force count_value to reach max

      // =====================================================
      // 4. RANDOM STRESS PHASE
      // =====================================================
      repeat (500) begin
        tr_gen = new();
        assert(tr_gen.randomize());
        gen_to_dr.put(tr_gen);
      end

      // Idle cycles
      repeat (5)
        send(1, 0, 0, 5);

    endtask

  endclass
endpackage
