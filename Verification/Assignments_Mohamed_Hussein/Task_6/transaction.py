class Transaction:

    def __init__(self, a=0, b=0, op=0, rst_n=1, out = 0, c = 0):
        self.a = a
        self.b = b
        self.op = op
        self.rst_n = rst_n
        self.out = out
        self.c = c

    def display(self, prefix=""):
        print(f"{prefix} a={self.a} b={self.b} op={self.op} rst_n={self.rst_n} out={self.out}")