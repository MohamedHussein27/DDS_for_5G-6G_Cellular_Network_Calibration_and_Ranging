package alu_scoreboard_pkg;
  
import alu_pkg::*;
import alu_sequnce_item_pkg::*;

class alu_scoreboard;
  

  mailbox #(alu_seq_item) mon2sb_mb;

  function new();
    
  endfunction
  
  

  // need to make error counter in real scoreboard
  int error_count = 0;


  // Reference model 
  /*
  function bit [4:0] calc_expected(
    bit [3:0] a,
    bit [3:0] b,
    opcode_t  op
  );
  
    case (op)
      ADD_op: return a + b;
      XOR_op: return a ^ b;
      AND_op: return a & b;
      OR_op : return a | b;
      default: return 'x;
    endcase
   // end
  endfunction
*/
  task run();
    alu_seq_item item;
    bit [4:0] expected;

    forever begin
      mon2sb_mb.get(item);

    //  expected = calc_expected(item.a, item.b, item.op);
        case (item.op)
      ADD_op:expected=(item.a +  item.b);
      XOR_op:expected=(item.a ^  item.b);
      AND_op:expected =(item.a & item.b);
      OR_op :expected= (item.a | item.b);
      
    endcase


      if (item.result !== expected) begin
         error_count++;
        $display(" ERROR A=%0d B=%0d OP=%0d EXP=%0d GOT=%0d errors=%0d",
                  item.a, item.b, item.op,
                  expected, item.result, error_count);
                 
      end
      else
        $display(" PASS  A=%0d B=%0d OP=%0d RESULT=%0d",
                  item.a, item.b, item.op, item.result);
    end
  endtask

endclass
endpackage