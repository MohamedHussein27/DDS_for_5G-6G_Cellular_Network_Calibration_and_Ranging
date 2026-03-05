class ALUTransaction:
    def __init__(self, a, b, op, rst_n=1):
        self.rst_n = rst_n
        self.a = a
        self.b = b
        self.op = op
        self.c = None
        self.out = None

    def __str__(self):
        return (f"ALUTransaction(rst_n={self.rst_n}, "
                f"ALUTransaction(a={self.a}, "
                f"b={self.b}, op={self.op}, "
                f"out={self.out}, c={self.c})")
