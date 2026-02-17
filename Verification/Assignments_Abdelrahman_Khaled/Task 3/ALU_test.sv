package ALU_test_pkg;

    import uvm_pkg::*;
    `include "uvm_macros.svh"

    import ALU_env_pkg::*;
    import ALU_config_pkg::*;

    import ALU_seq_item_pkg::*;
    import ALU_add_sequence_pkg::*;
    import ALU_xor_sequence_pkg::*;
    import ALU_and_sequence_pkg::*;
    import ALU_or_sequence_pkg::*;
    import ALU_random_sequence_pkg::*;

    //----------------------------------------------------------------------------
    // Base ALU Test
    //----------------------------------------------------------------------------
    // Responsibilities:
    // - Create environment
    // - Create and configure ALU_config
    // - Provide virtual interface to components
    //----------------------------------------------------------------------------
    class ALU_test extends uvm_test;

    `uvm_component_utils(ALU_test)

    // Environment handle
    ALU_env env;

    // Configuration object
    ALU_config config_test;

    // Virtual interface
    virtual ALU_if vif;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_test", uvm_component parent = null);
        super.new(name, parent);
    endfunction

    //--------------------------------------------------------------------------
    // Build Phase
    //--------------------------------------------------------------------------
    function void build_phase(uvm_phase phase);
        super.build_phase(phase);

        // Create environment
        env = ALU_env::type_id::create("env", this);

        // Create configuration object
        config_test = ALU_config::type_id::create("config_test");

        // Get virtual interface from top-level TB
        if (!uvm_config_db#(virtual ALU_if)::get(this, "", "ALU_IF", config_test.vif))
        `uvm_fatal("ALU_TEST", "Unable to get ALU virtual interface")

        // Set config object for all components below
        uvm_config_db#(ALU_config)::set(this, "*", "CFG", config_test);
    endfunction


    //--------------------------------------------------------------------------
    // End of elaboration Phase (show hierarchy)
    //--------------------------------------------------------------------------
    function void end_of_elaboration_phase(uvm_phase phase);
        super.end_of_elaboration_phase(phase);
        
        // This prints the entire UVM hierarchy in a tree format
        uvm_top.print_topology(); 
        
      endfunction

    endclass : ALU_test

    //----------------------------------------------------------------------------
    // ADD Test
    //----------------------------------------------------------------------------
    // Runs only ADD directed sequence
    //----------------------------------------------------------------------------
    class add_test extends ALU_test;

    `uvm_component_utils(add_test)

    ALU_add_sequence seq;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "add_test", uvm_component parent = null);
        super.new(name, parent);
    endfunction

    //--------------------------------------------------------------------------
    // Run Phase
    //--------------------------------------------------------------------------
    task run_phase(uvm_phase phase);
        super.run_phase(phase);

        // Create sequence
        seq = ALU_add_sequence::type_id::create("seq");

        // Start test
        phase.raise_objection(this);
        `uvm_info("ADD_TEST", "ADD sequence started", UVM_LOW)

        // Start sequence on ALU sequencer
        seq.start(env.ag_env.sr_ag);

        `uvm_info("ADD_TEST", "ADD sequence finished", UVM_LOW)
        phase.drop_objection(this);
    endtask

    endclass : add_test

    //----------------------------------------------------------------------------
    // XOR Test
    //----------------------------------------------------------------------------
    class xor_test extends ALU_test;

    `uvm_component_utils(xor_test)

    ALU_xor_sequence seq;

    function new(string name = "xor_test", uvm_component parent = null);
        super.new(name, parent);
    endfunction

    task run_phase(uvm_phase phase);
        super.run_phase(phase);

        seq = ALU_xor_sequence::type_id::create("seq");

        phase.raise_objection(this);
        `uvm_info("XOR_TEST", "XOR sequence started", UVM_LOW)

        seq.start(env.ag_env.sr_ag);

        `uvm_info("XOR_TEST", "XOR sequence finished", UVM_LOW)
        phase.drop_objection(this);
    endtask

    endclass

    //----------------------------------------------------------------------------
    // AND Test
    //----------------------------------------------------------------------------
    class and_test extends ALU_test;

    `uvm_component_utils(and_test)

    ALU_and_sequence seq;

    function new(string name = "and_test", uvm_component parent = null);
        super.new(name, parent);
    endfunction

    task run_phase(uvm_phase phase);
        super.run_phase(phase);

        seq = ALU_and_sequence::type_id::create("seq");

        phase.raise_objection(this);
        `uvm_info("AND_TEST", "AND sequence started", UVM_LOW)

        seq.start(env.ag_env.sr_ag);

        `uvm_info("AND_TEST", "AND sequence finished", UVM_LOW)
        phase.drop_objection(this);
    endtask

    endclass

    //----------------------------------------------------------------------------
    // OR Test
    //----------------------------------------------------------------------------
    class or_test extends ALU_test;

    `uvm_component_utils(or_test)

    ALU_or_sequence seq;

    function new(string name = "or_test", uvm_component parent = null);
        super.new(name, parent);
    endfunction

    task run_phase(uvm_phase phase);
        super.run_phase(phase);

        seq = ALU_or_sequence::type_id::create("seq");

        phase.raise_objection(this);
        `uvm_info("OR_TEST", "OR sequence started", UVM_LOW)

        seq.start(env.ag_env.sr_ag);

        `uvm_info("OR_TEST", "OR sequence finished", UVM_LOW)
        phase.drop_objection(this);
    endtask

    endclass

    //----------------------------------------------------------------------------
    // Random Test
    //----------------------------------------------------------------------------
    class random_test extends ALU_test;

    `uvm_component_utils(random_test)

    ALU_random_sequence seq;

    function new(string name = "random_test", uvm_component parent = null);
        super.new(name, parent);
    endfunction

    task run_phase(uvm_phase phase);
        super.run_phase(phase);

        seq = ALU_random_sequence::type_id::create("seq");

        phase.raise_objection(this);
        `uvm_info("RANDOM_TEST", "Random sequence started", UVM_LOW)

        seq.start(env.ag_env.sr_ag);

        `uvm_info("RANDOM_TEST", "Random sequence finished", UVM_LOW)
        phase.drop_objection(this);
    endtask

    endclass

endpackage : ALU_test_pkg
