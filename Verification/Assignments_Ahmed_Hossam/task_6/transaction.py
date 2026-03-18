import cocotb 
from cocotb_coverage.crv import *
class  transaction (Randomized):
    def __init__(self ,name = "TRANSACTION"):
        Randomized.__init__(self)
        self.name = name
        self.reset         =  0 
        self.a           =  0
        self.b           =  0
        self.op     =  0
        self.c =  0 
        self.out =  0 
        
        self.add_rand("reset"        , list(range(0,2)       )   ) 
        self.add_rand("a"          , list(range(0,16)   )   ) 
        self.add_rand("b"          , list(range(0,16)   )   )
        self.add_rand("op"    , list(range(0,4)      )   )


        
    def display(self,name = "TRANSACTION"):
        cocotb.log.info("******************"+str(name)+"*******************")
        cocotb.log.info("the Value of reset        is   " + str(self.reset        ))
        cocotb.log.info("the Value of a           is   "  + str(self.a         ))
        cocotb.log.info("the Value of b           is   " + str(self.b          ))
        cocotb.log.info("the Value of op     is   " + str(self.op    ))
        cocotb.log.info("the Value of c is   " +   str(self.c ))
        cocotb.log.info("the Value of out  is   " + str(self.out ))
        cocotb.log.info("**************************************************")
