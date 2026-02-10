import shared_pkg::*;
package counter_seq_item_pkg;
    class counter_seq_item;
        bit clk;
        rand bit rst_n;
        rand bit start;
        rand bit [15:0] wait_timer;
        rand bit flag;
        logic busy;
        logic [4:0] count_value;

        // constraints
        constraint rst_n_con {
            rst_n dist {0 := 1, 1 := 200}; // reset active low
        }

        constraint flag_con {
            flag dist {0 := 95, 1 := 5}; // higher values has less prop.
        }

        constraint start_con {
            start dist {0 := 90, 1 := 10}; // higher values has less prop.
        }
        
    endclass
endpackage