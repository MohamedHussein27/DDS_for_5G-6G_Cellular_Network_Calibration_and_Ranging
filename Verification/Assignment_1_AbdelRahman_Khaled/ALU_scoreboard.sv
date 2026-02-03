package ALU_scoreboard_pkg;

  import ALU_transaction_pkg::*;

  // ALU Scoreboard class
  // Compares DUT output with golden reference model
  class ALU_scoreboard;

    // Reference signals
    logic       c_ref = 0;
    logic [3:0] out_ref = 0;

    // Enumeration for ALU operations
    typedef enum logic [1:0] {
      ADD_op = 0,
      XOR_op = 1,
      AND_op = 2,
      OR_op  = 3
    } operation_t;

    // Counters
    int Error_count   = 0;
    int Correct_count = 0;

    // Transaction from monitor
    ALU_transaction tr_sb;

    // Mailbox connected to monitor
    mailbox #(ALU_transaction) mbx_mon_to_sb;

    // Constructor
    function new(mailbox #(ALU_transaction) mbx);
      this.mbx_mon_to_sb = mbx;
    endfunction

    // Reference model
    function void ALU_ref(ALU_transaction tr_ref);
      unique case(tr_ref.op)
        ADD_op: {c_ref, out_ref} = tr_ref.a + tr_ref.b;
        XOR_op: begin out_ref = tr_ref.a ^ tr_ref.b; c_ref = 0; end
        AND_op: begin out_ref = tr_ref.a & tr_ref.b; c_ref = 0; end
        OR_op : begin out_ref = tr_ref.a | tr_ref.b; c_ref = 0; end
        default: begin c_ref = 0; out_ref = '0; end
      endcase
    endfunction

    // Check DUT transaction
    function void ALU_check(ALU_transaction tr_DUT);

      ALU_ref(tr_DUT);

      if ((tr_DUT.op == ADD_op && {tr_DUT.c, tr_DUT.out} !== {c_ref, out_ref}) ||
          (tr_DUT.op != ADD_op && (tr_DUT.out !== out_ref || tr_DUT.c !== 0))) begin
        Error_count++;
        $display("ERROR at time %0t: a=%0d b=%0d op=%0b DUT={c=%0b out=%0d} Expected={c=%0b out=%0d}",
                  $time, tr_DUT.a, tr_DUT.b, tr_DUT.op, tr_DUT.c, tr_DUT.out, c_ref, out_ref);
      end else begin
        Correct_count++;
      end

    endfunction : ALU_check

    // Main run task: consume from monitor mailbox
    task run();
      forever begin
        mbx_mon_to_sb.get(tr_sb);
        ALU_check(tr_sb);
      end
    endtask

    // Report at end of simulation
    task report();
      $display("Simulation summary: Correct=%0d Errors=%0d", Correct_count, Error_count);
    endtask

  endclass : ALU_scoreboard

endpackage : ALU_scoreboard_pkg
