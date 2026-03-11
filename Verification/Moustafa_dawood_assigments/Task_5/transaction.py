

# we first dont have interface in python
# we need only data container 

class ALUTransaction:

    def __init__(self, a=0, b=0, op=0):
        self.a = a
        self.b = b
        self.op = op
        self.result = None