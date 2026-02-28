class ALUTransaction:
    def __init__(self, a, b, op):
        self.a = a
        self.b = b
        self.op = op
        self.c = None
        self.out = None

    def __str__(self):
        return (f"ALUTransaction(a={self.a}, "
                f"b={self.b}, op={self.op}, "
                f"out={self.out}, c={self.c})")
