

# we first dont have interface in python
# we need only data container 

class ALUTransaction:

    def __init__(self, a, b, op):
        self.a = a
        self.b = b
        self.op = op
        self.result = None