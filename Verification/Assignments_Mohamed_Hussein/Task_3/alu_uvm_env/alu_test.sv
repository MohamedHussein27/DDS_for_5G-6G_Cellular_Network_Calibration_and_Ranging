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

    class alu_base_test extends uvm_test;
        `uvm_component_utils(alu_base_test)

        alu_env            env;
        alu_config_obj     alu_cfg;
        alu_reset_sequence reset_seq;

        function new(string name = "alu_base_test", uvm_component parent = null);
            super.new(name, parent);
        endfunction

        function void build_phase(uvm_phase phase);
            super.build_phase(phase);

            env       = alu_env::type_id::create("env", this);
            alu_cfg   = alu_config_obj::type_id::create("alu_cfg");
            reset_seq = alu_reset_sequence::type_id::create("reset_seq");

            if (!uvm_config_db #(virtual alu_if)::get(this,"","alu_V", alu_cfg.alu_vif))
            `uvm_fatal("build_phase", "Main interface not found");

            uvm_config_db #(alu_config_obj)::set(this,"*","CFG", alu_cfg);
        endfunction


        // (runs automatically before main_phase)
        task reset_phase(uvm_phase phase);
            super.reset_phase(phase);

            phase.raise_objection(this);

            `uvm_info("RESET_PHASE", "Starting Reset Sequence", UVM_LOW)
            reset_seq.start(env.agt.sqr);
            `uvm_info("RESET_PHASE", "Reset Completed", UVM_LOW)

            // debug when it excutes
            `uvm_info("RESET_PHASE", $sformatf("Time = %0t", $time), UVM_LOW) 

            phase.drop_objection(this);
        endtask

    endclass



    class alu_add_xor_test extends alu_base_test;
        `uvm_component_utils(alu_add_xor_test)

        alu_add_sequence add_seq;
        alu_xor_sequence xor_seq;

        function new(string name = "alu_add_xor_test", uvm_component parent = null);
            super.new(name, parent);
        endfunction

        function void build_phase(uvm_phase phase);
            super.build_phase(phase);

            add_seq = alu_add_sequence::type_id::create("add_seq");
            xor_seq = alu_xor_sequence::type_id::create("xor_seq");
        endfunction


        task main_phase(uvm_phase phase);
            phase.raise_objection(this);

            add_seq.start(env.agt.sqr);
            xor_seq.start(env.agt.sqr);

            phase.drop_objection(this);
        endtask

    endclass


    class alu_and_or_test extends alu_base_test;
        `uvm_component_utils(alu_and_or_test)    
        
        
        // Specific Directed Sequences based on your imports
        alu_and_sequence   and_seq;
        alu_or_sequence    or_seq;

        function new(string name = "alu_and_or_test", uvm_component parent = null);
            super.new(name, parent);
        endfunction

        function void build_phase(uvm_phase phase);
            super.build_phase(phase);

            // Sequence Creation
            
            and_seq   = alu_and_sequence::type_id::create("and_seq");
            or_seq    = alu_or_sequence::type_id::create("or_seq");
        endfunction

        task main_phase(uvm_phase phase);
            super.run_phase(phase);
            phase.raise_objection(this);

            `uvm_info("run_phase", "Starting AND Sequence", UVM_LOW)
            and_seq.start(env.agt.sqr);

            `uvm_info("run_phase", "Starting OR Sequence", UVM_LOW)
            or_seq.start(env.agt.sqr);

            phase.drop_objection(this);
        endtask
    endclass

    class alu_rand_test extends alu_base_test;
        `uvm_component_utils(alu_rand_test)

        alu_main_sequence  main_seq;

        function new(string name = "alu_rand_test", uvm_component parent = null);
            super.new(name, parent);
        endfunction

        function void build_phase(uvm_phase phase);
            super.build_phase(phase);

            main_seq  = alu_main_sequence::type_id::create("main_seq");
        endfunction

        task main_phase(uvm_phase phase);
            phase.raise_objection(this);

            main_seq.start(env.agt.sqr);

            phase.drop_objection(this);
        endtask
    endclass

endpackage