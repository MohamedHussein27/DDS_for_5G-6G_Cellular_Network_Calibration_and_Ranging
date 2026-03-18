from transaction import transaction
from cocotb.triggers import *
import cocotb
import cocotb.queue
from cocotb_coverage.coverage import CoverPoint, coverage_db

@CoverPoint("top.out", vname="out", bins=list(range(0, 16))) 
@CoverPoint("top.c",   vname="c",   bins=list(range(0, 2)))
@CoverPoint("top.op",  vname="op",  bins=list(range(0, 4)))
def sample(out, c, op):
    pass

class subscriber():
   
    def __init__(self, name="SUBSCRIBER"): 
        self.name     = name
    
        self.sub_mail = cocotb.queue.Queue()

    async def run_subscriber(self): 
        cocotb.log.info("[Subscriber] STARTING.")
        
        while True:
            # Wait for the monitor to drop a packet in the queue
            self.t_sub = await self.sub_mail.get() 
            
            cocotb.log.debug("[Subscriber] Packet received from monitor. Sampling...") 
            
            # Pass the individual variables into the sample function
            sample(
                self.t_sub.out, 
                self.t_sub.c, 
                self.t_sub.op
            )
 
    def coverage_report(self):
        # Fetch the data using the exact names from the @CoverPoint decorators
        out_cov = coverage_db["top.out"].coverage          
        out_p   = coverage_db["top.out"].cover_percentage  
        
        cocotb.log.info(f"The 'out' coverage is : {out_cov}")
        cocotb.log.info(f"The 'out' coverage percentage is : {out_p}%")
        
        coverage_db.export_to_xml("alu_coverage.xml")