package ALU_scoreboard_pkg;

  // Import transaction definition
  import ALU_transaction_pkg::*;

  // Import shared package definition
  import ALU_shared_pkg::*;

  // Reference signals
  logic       c_ref;
  logic [3:0] out_ref;

  // Enumeration for ALU operations
  typedef enum logic [1:0] {
    ADD_op = 0,
    XOR_op = 1,
    AND_op = 2,
    OR_op  = 3
  } operation_t;

  // ALU Scoreboard class
  class ALU_scoreboard;

    // Reference model (golden model)
    // Computes the expected output based on inputs
    function void ALU_ref(ALU_transaction tr_ref);
      case (tr_ref.op)
        ADD_op: {c_ref, out_ref} = tr_ref.a + tr_ref.b;
        XOR_op: out_ref         = tr_ref.a ^ tr_ref.b;
        AND_op: out_ref         = tr_ref.a & tr_ref.b;
        OR_op : out_ref         = tr_ref.a | tr_ref.b;
        default: begin
          c_ref   = 0;
          out_ref = '0;
        end
      endcase
    endfunction : ALU_ref

    // Check DUT output against expected result
    function void ALU_check(ALU_transaction tr_DUT);

        // Compute golden/reference values
        ALU_ref(tr_DUT);

        // Compare DUT output with golden model
        if ((tr_DUT.op == ADD_op) && ({tr_DUT.c, tr_DUT.out} != {c_ref, out_ref})) begin
            Error_count++;
            $display("ERROR: ADD operation failed at time %0t", $time);
            $display("       a=%0d, b=%0d, DUT={c=%0b, out=%0d}, Expected={c=%0b, out=%0d}",
                    tr_DUT.a, tr_DUT.b, tr_DUT.c, tr_DUT.out, c_ref, out_ref);

        end else if ((tr_DUT.op == XOR_op) && (tr_DUT.out != out_ref)) begin
            Error_count++;
            $display("ERROR: XOR operation failed at time %0t", $time);
            $display("       a=%0d, b=%0d, DUT out=%0d, Expected out=%0d",
                    tr_DUT.a, tr_DUT.b, tr_DUT.out, out_ref);

        end else if ((tr_DUT.op == AND_op) && (tr_DUT.out != out_ref)) begin
            Error_count++;
            $display("ERROR: AND operation failed at time %0t", $time);
            $display("       a=%0d, b=%0d, DUT out=%0d, Expected out=%0d",
                    tr_DUT.a, tr_DUT.b, tr_DUT.out, out_ref);

        end else if ((tr_DUT.op == OR_op) && (tr_DUT.out != out_ref)) begin
            Error_count++;
            $display("ERROR: OR operation failed at time %0t", $time);
            $display("       a=%0d, b=%0d, DUT out=%0d, Expected out=%0d",
                    tr_DUT.a, tr_DUT.b, tr_DUT.out, out_ref);

        end else begin
            Correct_count++;
        end

    endfunction : ALU_check

  endclass : ALU_scoreboard

endpackage : ALU_scoreboard_pkg

  