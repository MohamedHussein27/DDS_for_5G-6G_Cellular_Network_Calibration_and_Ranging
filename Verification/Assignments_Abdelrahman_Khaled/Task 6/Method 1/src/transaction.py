import random


class ALUTransaction:
    def __init__(self, a=0, b=0, op=0, rst_n=1):
        self.rst_n = rst_n
        self.a = a
        self.b = b
        self.op = op
        self.c = None
        self.out = None

    def randomize(self):
        """Randomize all inputs."""
        self.a = random.randint(0, 15)
        self.b = random.randint(0, 15)
        self.op = random.randint(0, 3)
        self.rst_n = 1

    def randomize_with(self, constraint):
        """
        Randomize until constraint(self) returns True.
        Example:
            txn.randomize_with(lambda t: t.op == 0 and t.a < t.b)
        """
        for _ in range(1000):
            self.randomize()
            if constraint(self):
                return
        raise ValueError("Constraint could not be satisfied")

    def __str__(self):
        return (f"ALUTransaction(rst_n={self.rst_n}, "
                f"a={self.a}, b={self.b}, op={self.op}, "
                f"out={self.out}, c={self.c})")
