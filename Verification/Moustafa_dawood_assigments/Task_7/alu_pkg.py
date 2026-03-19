from enum import IntEnum

class opcode_t(IntEnum):
    ADD_op = 0
    XOR_op = 1
    AND_op = 2
    OR_op  = 3

class error_packet_t:
    def __init__(self, a=0, b=0, op=0):
        self.a = a
        self.b = b
        self.op = op