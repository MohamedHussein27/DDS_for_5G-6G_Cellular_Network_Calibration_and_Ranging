package ALU_generator_pkg;

  // Import ALU transaction definition
  import ALU_transaction_pkg::*;

  // ALU generator class
  // Generates exhaustive stimulus for all ALU operations
  class ALU_generator;

    // Transaction object used for stimulus generation
    ALU_transaction tr_gen;

    // Mailbox to send transactions to the driver
    mailbox #(ALU_transaction) gen_to_dr;

    // Constructor
    function new(mailbox #(ALU_transaction) gen_to_dr);
      this.gen_to_dr = gen_to_dr;
      tr_gen = new();
    endfunction : new

    // Main stimulus generation task
    task run();

      // Loop over all ALU operations
      for (int op = 0; op < 4; op++) begin

        // Loop over all possible values of operand a
        for (int a = 0; a < 16; a++) begin

          // Loop over all possible values of operand b
          for (int b = 0; b < 16; b++) begin

            // Create a new transaction for each combination
            tr_gen = new();

            // Assign transaction fields
            tr_gen.op = op[1:0]; // cast int to 2-bit opcode
            tr_gen.a  = a[3:0];
            tr_gen.b  = b[3:0];

            // Send transaction to driver
            gen_to_dr.put(tr_gen);

          end
        end
      end

    endtask : run

  endclass : ALU_generator

endpackage : ALU_generator_pkg
