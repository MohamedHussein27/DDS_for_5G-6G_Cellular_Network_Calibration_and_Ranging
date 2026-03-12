from cocotb.triggers import *
from transaction import *
import cocotb
import cocotb.queue

class monitor ():
   t_monitor    = transaction()
   mon_mail_s   = cocotb.queue.Queue()
   mon_mail_su  = cocotb.queue.Queue()
   def __int__(self,name = "MONITOR"):
        self.name= name
      
   async def run_monitor (self,dut_monitor):
      cocotb.log.info("[Monitor] STARTING.")
      await RisingEdge(dut_monitor.clk)
      while(True):
        cocotb.log.info("[Monitor] waiting for item ...")
        await RisingEdge(dut_monitor.clk)
        await ReadOnly()
        self.t_monitor.reset            =   int(dut_monitor.reset.value )
        self.t_monitor.a              =   int(dut_monitor.a.value )    
        self.t_monitor.b              =   int(dut_monitor.b.value )    
        self.t_monitor.op        =   int(dut_monitor.op.value  )
        self.t_monitor.c     =   int(dut_monitor.c.value )
        self.t_monitor.out     =   int(dut_monitor.out.value  )
        await self.mon_mail_s.put(self.t_monitor)  
        await self.mon_mail_su.put(self.t_monitor)