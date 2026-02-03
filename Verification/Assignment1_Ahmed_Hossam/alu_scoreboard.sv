`include "alu_item.sv"
class alu_scoreboard;
    mailbox #(alu_item) mon2scb;
    int error_count = 0;
    int pass_count = 0;

    function new(mailbox mon2scb);
        this.mon2scb = mon2scb;
    endfunction

    task run();
        forever begin
            alu_item item;
            mon2scb.get(item);
            check(item);
        end
    endtask

    task check(alu_item item);
        bit [4:0] expected_result;
        bit expected_c;
        bit [3:0] expected_out;

        // Golden Model Logic
        case (item.op)
            ADD_op: expected_result = item.a + item.b;
            XOR_op: expected_result = {1'b0, item.a ^ item.b};
            AND_op: expected_result = {1'b0, item.a & item.b};
            OR_op:  expected_result = {1'b0, item.a | item.b};
            default: expected_result = 0;
        endcase

        expected_c   = expected_result[4];
        expected_out = expected_result[3:0];

        if (item.c != expected_c || item.out != expected_out) begin
            error_count++;
            $display("ERROR A=%0d B=%0d op=%0d and expected_c =%0d and expected_out =%0d and c = %0d and out =%0d",
                  item.a, item.b, item.op,
                  expected_c, expected_out , item.c , item.out , $time );
                  $display("ERROR count =%0d", error_count);
        end else begin
            $display("PASS  A=%0d B=%0d OP=%0d expected_c =%0d and expected_out =%0d and c = %0d and out =%0d",
                  item.a, item.b, item.op,expected_c, expected_out , item.c , item.out , $time );
            pass_count++;
             $display("PASS count =%0d", pass_count);

        end
    endtask

    
endclass