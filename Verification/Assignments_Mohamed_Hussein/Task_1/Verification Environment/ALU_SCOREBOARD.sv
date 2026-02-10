import shared_pkg::*;
import alu_seq_item_pkg::*;
class alu_scoreboard;


    // reference signals
    logic [3:0] out_ref;
    logic c_ref;

    mailbox #(alu_seq_item) mon2sb;

    function new();
    endfunction

    function void connect_mail(mailbox #(alu_seq_item) m);
        mon2sb = m;
    endfunction

    // compare function
    function void check_data (alu_seq_item tr);
        reference_model(tr);
        // compare
        if (tr.out !== out_ref) begin
            error_count_out++;
            $display("error in data out, ref_out is: %0d     while dut out is: %0d", out_ref, tr.out);
        end
        else
            correct_count_out++;
        if(tr.c !== c_ref) begin
            error_count_c++;
            $display("error in carry, ref_c is: %0d    while dut c is: %0d", c_ref, tr.c);
        end
        else
            correct_count_c++;
    endfunction

    // reference function
    function void reference_model (alu_seq_item tr_ref);
        if (~tr_ref.rst_n) 
            {c_ref, out_ref} = 5'b0;
        else
            case (tr_ref.op)
                2'b00: begin // addition
                    {c_ref, out_ref} = tr_ref.a + tr_ref.b;
                end
                2'b01: begin // XOR
                    c_ref = 0;
                    out_ref = tr_ref.a ^ tr_ref.b;
                end
                2'b10: begin // AND
                    c_ref = 0;
                    out_ref = tr_ref.a &  tr_ref.b;
                end
                2'b11: begin // OR
                    c_ref = 0;
                    out_ref = tr_ref.a | tr_ref.b;
                end
                default: begin
                    c_ref = 0;
                    out_ref = 4'b0000;
                end
            endcase
    endfunction

    task run();
        alu_seq_item item;
        forever begin
            mon2sb.get(item);
            check_data(item);
        end
    endtask

endclass
