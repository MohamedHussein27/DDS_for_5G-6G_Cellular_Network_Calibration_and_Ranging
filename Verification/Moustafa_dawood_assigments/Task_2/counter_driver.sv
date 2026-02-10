/*
package counter_driver_pkg;
import counter_seq_item_pkg::*;

class counter_driver;

  virtual counter_bfm bfm;
  mailbox #(counter_seq_item) seq2drv_mb;

  function new();
  endfunction

  // --------------------------------------
  // Reset task
  // --------------------------------------
  task drive_reset();
    $display("DRIVER Asserting Reset");
    bfm.rst_n <= 0;                 // Drive reset low
    repeat (2) @(posedge bfm.clk);  // Hold for 2 clocks
    bfm.rst_n <= 1;                 // Release reset
    @(posedge bfm.clk);
    $display("DRIVER Reset Done");
  endtask

  // --------------------------------------
  // Main driver loop (drive at negedge clk)
  // --------------------------------------
  task run();
    counter_seq_item item;
    forever begin
      seq2drv_mb.get(item);

      // Drive signals at negedge clk
      @(negedge bfm.clk);
      bfm.start      <= item.start;
      bfm.flag       <= item.flag;
      bfm.rst_n      <= item.rst_n;
      bfm.wait_timer <= item.wait_timer;

      // Hold them for one full clock to ensure DUT samples correctly
      @(negedge bfm.clk);
    end
  endtask

endclass

endpackage
*/
package counter_driver_pkg;
  import counter_seq_item_pkg::*;

  class counter_driver;
    virtual counter_bfm bfm;
    mailbox #(counter_seq_item) seq2drv_mb;

    function new();
    endfunction

    // ---------------------------------
    // Reset
    // ---------------------------------
    task drive_reset();
      $display("DRIVER Asserting Reset");
      bfm.rst_n      <= 0;
      bfm.start      <= 0;
      bfm.flag       <= 0;
      bfm.wait_timer <= 0;

      repeat (2) @(posedge bfm.clk);
      bfm.rst_n <= 1;
      @(posedge bfm.clk);
      $display("DRIVER Reset Done");
    endtask

    // ---------------------------------
    // Driver main loop
    // ---------------------------------
    task run();
      counter_seq_item item;

      forever begin
        seq2drv_mb.get(item);

        @(negedge bfm.clk);
        bfm.rst_n      <= item.rst_n;
        bfm.start      <= item.start;
        bfm.flag       <= item.flag;
        bfm.wait_timer <= item.wait_timer;
      end
    endtask

  endclass
endpackage
