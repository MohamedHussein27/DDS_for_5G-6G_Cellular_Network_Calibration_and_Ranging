from transaction    import * 
from generator       import * 
from monitor         import *
from driver     import *
from scoreboard import *
from subscriber import subscriber 
import cocotb
import cocotb.queue


class environment() :
   
    def __init__(self ,name = "ENVIRONMENT"): 
        self.name        = name
        self.join_any    = Event(name=None) 
        self.g           = generator(self.join_any )
        self.d           = driver()
        self.m           = monitor()
        self.s           = scoreboard()
        self.su          = subscriber()                  

    async def run_environment (self,dut) : 
        cocotb.log.info(" WE ARE STARTING THE SIMULATION ")
        self.d.driv_mail     = self.g.gen_mail
        self.d.driv_handover = self.g.gen_handover
        self.s.score_mail    = self.m.mon_mail_s
        self.su.sub_mail     = self.m.mon_mail_su 
        cocotb.start_soon(self.g.run_generator(dut))
        cocotb.start_soon(self.d.run_driver(dut)   )
        cocotb.start_soon(self.m.run_monitor(dut)  )
        cocotb.start_soon(self.s.run_scoreboard()  )
        cocotb.start_soon(self.su.run_subscriber() )