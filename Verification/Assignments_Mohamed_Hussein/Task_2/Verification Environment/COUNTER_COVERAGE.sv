import counter_seq_item_pkg::*;
class counter_coverage;
    string name;

    mailbox #(counter_seq_item) mon2cov;

    

    counter_seq_item cov_item;

    // cover group
    covergroup counter_Cross_Group;
        // cover points
        START_CP: coverpoint cov_item.start;

        FLAG_CP: coverpoint cov_item.flag;

        BUSY_CP: coverpoint cov_item.busy;

        COUNT_CP: coverpoint cov_item.count_value;

    endgroup
    
    function new(string name = "SUBSCRIBER");
        this.name = name;
        mon2cov = new();
        counter_Cross_Group = new();
    endfunction

    task run();
        counter_seq_item item;
        
        forever begin
            item = new();
            mon2cov.get(item);
            cov_item = item;
            counter_Cross_Group.sample();
        end
    endtask

endclass