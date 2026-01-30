import alu_driver_pkg::*;
import alu_monitor_pkg::*;
import alu_seq_item_pkg::*;
import shared_pkg::*;
module alu_top;

    alu_if aluif();

    ALU dut (aluif);

    alu_driver  drv;
    alu_monitor mon;
    alu_seq_item tr;

    initial begin
        drv = new(aluif);
        mon = new(aluif);
        tr  = new();

        repeat (10000) begin
            assert(tr.randomize());
            drv.drive(tr);
            mon.sample(tr);
        end
        $display("error count = %0d, correct count = %0d", error_count_out, correct_count_out);
        $stop;
    end

endmodule