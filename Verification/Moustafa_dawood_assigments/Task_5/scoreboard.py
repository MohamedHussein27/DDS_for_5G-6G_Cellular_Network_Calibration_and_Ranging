import logging


class ALUScoreboard:

    def __init__(self):

        self.mon2sb = None
        self.error_count = 0

    async def run(self):

        while True:

            item = await self.mon2sb.get()

            if item.op == 0:
                expected = item.a + item.b
            elif item.op == 1:
                expected = item.a ^ item.b
            elif item.op == 2:
                expected = item.a & item.b
            elif item.op == 3:
                expected = item.a | item.b

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