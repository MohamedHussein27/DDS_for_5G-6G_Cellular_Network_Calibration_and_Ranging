# scoreboard.py
import logging


class ALUScoreboard:

    def __init__(self):
        self.error_count = 0

    def check(self, item):

        # Reference model
        if item.op == 0:          # ADD_op
            expected = item.a + item.b

        elif item.op == 1:        # XOR_op
            expected = item.a ^ item.b

        elif item.op == 2:        # AND_op
            expected = item.a & item.b

        elif item.op == 3:        # OR_op
            expected = item.a | item.b

        else:
            raise ValueError("Invalid opcode")

        if item.result != expected:

            self.error_count += 1

            logging.error(
                f"ERROR A={item.a} "
                f"B={item.b} "
                f"OP={item.op} "
                f"EXP={expected} "
                f"GOT={item.result} "
                f"errors={self.error_count}"
            )

        else:

            logging.info(
                f"PASS A={item.a} "
                f"B={item.b} "
                f"OP={item.op} "
                f"RESULT={item.result}"
            )