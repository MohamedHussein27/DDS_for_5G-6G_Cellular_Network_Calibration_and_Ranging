package alu_scoreboard_pkg;
    import alu_seq_item_pkg::*;
    import alu_shared_pkg::*;
    import uvm_pkg::*;
    `include "uvm_macros.svh"
    class alu_scoreboard extends uvm_scoreboard;
        `uvm_component_utils(alu_scoreboard)
        uvm_analysis_export #(alu_seq_item) sb_export;
        uvm_tlm_analysis_fifo #(alu_seq_item) sb_fifo;
        alu_seq_item seq_item_sb;


        // reference signals
        logic [3:0] out_ref;
        logic c_ref;

        function new(string name = "alu_scoreboard", uvm_component parent = null);
            super.new(name, parent);
        endfunction


        function void build_phase(uvm_phase phase);
            super.build_phase(phase);
            sb_export = new("sb_export", this);
            sb_fifo = new("sb_fifo", this);
        endfunction

        // connect
        function void connect_phase(uvm_phase phase);
            super.connect_phase(phase);
            sb_export.connect(sb_fifo.analysis_export);
        endfunction

        // run
        task run_phase(uvm_phase phase);
            super.run_phase(phase);
            forever begin
                sb_fifo.get(seq_item_sb);
                ref_model(seq_item_sb);

                if (seq_item_sb.out !== out_ref) begin
                    $display("error in out at %0d", error_count_out+1);
                    `uvm_error("run_phase", $sformatf("comparison out transaction received by the DUT:%s while the reference out:0b%0b",
                    seq_item_sb.convert2string(), out_ref));
                end
                else begin
                    `uvm_info("run_phase", $sformatf("correct out out: %s", seq_item_sb.convert2string()), UVM_HIGH);
                    correct_count_out++;
                end

                if (seq_item_sb.c !== c_ref) begin
                    $display("error in c at %0d", error_count_c+1);
                    `uvm_error("run_phase", $sformatf("comparison c transaction received by the DUT:%s while the reference out:0b%0b",
                    seq_item_sb.convert2string(), c_ref));
                end
                else begin
                    `uvm_info("run_phase", $sformatf("correct c : %s", seq_item_sb.convert2string()), UVM_HIGH);
                    correct_count_c++;
                end
            end
        endtask

        // reference model
        task ref_model (alu_seq_item seq_item_chk);
            if (~seq_item_chk.rst_n) 
                {c_ref, out_ref} = 5'b0;
            else
                case (seq_item_chk.op)
                    2'b00: begin // addition
                        {c_ref, out_ref} = seq_item_chk.a + seq_item_chk.b;
                    end
                    2'b01: begin // XOR
                        c_ref = 0;
                        out_ref = seq_item_chk.a ^ seq_item_chk.b;
                    end
                    2'b10: begin // AND
                        c_ref = 0;
                        out_ref = seq_item_chk.a &  seq_item_chk.b;
                    end
                    2'b11: begin // OR
                        c_ref = 0;
                        out_ref = seq_item_chk.a | seq_item_chk.b;
                    end
                    default: begin
                        c_ref = 0;
                        out_ref = 4'b0000;
                    end
                endcase
        endtask
        
        // report
        function void report_phase(uvm_phase phase);
            super.report_phase(phase);
            `uvm_info("report_phase", $sformatf("total successful transactions: %0d", correct_count_out), UVM_MEDIUM);
            `uvm_info("report_phase", $sformatf("total failed transactions: %0d", error_count_out), UVM_MEDIUM);
        endfunction
    endclass
endpackage