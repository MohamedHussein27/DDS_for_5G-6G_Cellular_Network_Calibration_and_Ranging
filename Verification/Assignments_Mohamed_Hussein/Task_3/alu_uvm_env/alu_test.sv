package alu_test_pkg;
    import alu_config_obj_pkg::*;
    import alu_env_pkg::*;
    import alu_reset_sequence_pkg::*;
    import alu_main_sequence_pkg::*;
    
    // Directed Sequence Packages
    import alu_add_sequence_pkg::*;
    import alu_xor_sequence_pkg::*;
    import alu_and_sequence_pkg::*;
    import alu_or_sequence_pkg::*;

    import uvm_pkg::*;
    `include "uvm_macros.svh"

    class alu_test extends uvm_test;
        `uvm_component_utils(alu_test)

        alu_env env; 
        alu_config_obj alu_cfg;    
        
        // Sequence Handles
        alu_reset_sequence reset_seq;
        alu_main_sequence  main_seq;
        
        // Specific Directed Sequences based on your imports
        alu_add_sequence   add_seq;
        alu_xor_sequence   xor_seq;
        alu_and_sequence   and_seq;
        alu_or_sequence    or_seq;

        function new(string name = "alu_test", uvm_component parent = null);
            super.new(name, parent);
        endfunction

        function void build_phase(uvm_phase phase);
            super.build_phase(phase);
            env     = alu_env::type_id::create("env", this);
            alu_cfg = alu_config_obj::type_id::create("alu_cfg");

            // Sequence Creation
            reset_seq = alu_reset_sequence::type_id::create("reset_seq");
            main_seq  = alu_main_sequence::type_id::create("main_seq");
            
            // Create specific directed sequences
            add_seq   = alu_add_sequence::type_id::create("add_seq");
            xor_seq   = alu_xor_sequence::type_id::create("xor_seq");
            and_seq   = alu_and_sequence::type_id::create("and_seq");
            or_seq    = alu_or_sequence::type_id::create("or_seq");

            // Config DB gets
            if (!uvm_config_db #(virtual alu_if)::get(this,"","alu_V", alu_cfg.alu_vif))
                `uvm_fatal("build_phase", "Main interface not found");

            uvm_config_db #(alu_config_obj)::set(this,"*","CFG", alu_cfg);
        endfunction

        task run_phase(uvm_phase phase);
            super.run_phase(phase);
            phase.raise_objection(this);

            // 1. Initial Reset to start the simulation
            `uvm_info("run_phase", "Initial reset asserted", UVM_LOW)
            reset_seq.start(env.agt.sqr);
            `uvm_info("run_phase", "Initial reset deasserted", UVM_LOW)

            // 2. Executing Directed Sequences First
            `uvm_info("run_phase", "Executing Directed Test Battery", UVM_LOW)

            `uvm_info("run_phase", "Starting ADD Sequence", UVM_LOW)
            add_seq.start(env.agt.sqr);

            `uvm_info("run_phase", "Starting XOR Sequence", UVM_LOW)
            xor_seq.start(env.agt.sqr);

            `uvm_info("run_phase", "Starting AND Sequence", UVM_LOW)
            and_seq.start(env.agt.sqr);

            `uvm_info("run_phase", "Starting OR Sequence", UVM_LOW)
            or_seq.start(env.agt.sqr);

            // 3. SECOND RESET: Prepare for Random Testing
            // We reset here to ensure the DUT is clean before the random sequence starts
            `uvm_info("run_phase", "Intermediate reset for Random Tests", UVM_LOW)
            reset_seq.start(env.agt.sqr);
            `uvm_info("run_phase", "Intermediate reset deasserted", UVM_LOW)

            // 4. Main Random/Standard Sequence (Run LAST)
            `uvm_info("run_phase", "Starting Main (Random) Sequence", UVM_LOW)
            main_seq.start(env.agt.sqr);
            `uvm_info("run_phase", "Main Sequence finished", UVM_LOW)

            phase.drop_objection(this);
        endtask
    endclass
endpackage