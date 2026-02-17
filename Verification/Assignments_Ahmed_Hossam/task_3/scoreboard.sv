`include "seq_item.sv"
class scoreboard extends uvm_scoreboard;
`uvm_component_utils(scoreboard)
uvm_analysis_export #(seq_item) score_export;
uvm_tlm_analysis_fifo #(seq_item) score_fifo;
seq_item item;
bit       c;
bit [3:0] out;
int error_count=0;
int correct_count=0;
function new(string name = "scoreboard", uvm_component parent = null);
super.new(name, parent);
endfunction
function void build_phase(uvm_phase phase);
super.build_phase(phase);
score_export = new("score_export", this);
score_fifo = new("score_fifo", this);
endfunction
function void connect_phase(uvm_phase phase);
super.connect_phase(phase);
score_export.connect(score_fifo.analysis_export);
endfunction
task run_phase(uvm_phase phase);
super.run_phase(phase);
forever begin 
    score_fifo.get(item);
    ref_model(item);
    if (item.out !== out || item.c !== c) begin
        `uvm_error("Scoreboard", $sformatf("Mismatch: Got out=%b, c=%b; Expected out=%b, c=%b", item.out, item.c, out, c));
        error_count++;
    end else begin
        `uvm_info("Scoreboard", $sformatf("Match: Got out=%b, c=%b; Expected out=%b, c=%b", item.out, item.c, out, c), UVM_HIGH);
        correct_count++;
    end
end
endtask
task ref_model(seq_item item);
case (item.op)
            ADD_op: {c,out} = item.a + item.b;
            XOR_op: {c,out} = {1'b0, item.a ^ item.b};
            AND_op: {c,out} = {1'b0, item.a & item.b};
            OR_op:  {c,out} = {1'b0, item.a | item.b};
            default: {c,out} = 0;
        endcase
endtask
function void report_phase(uvm_phase phase);
super.report_phase(phase);
`uvm_info("Scoreboard", $sformatf("Total Correct: %0d, Total Errors: %0d", correct_count, error_count), UVM_LOW);
endfunction
endclass