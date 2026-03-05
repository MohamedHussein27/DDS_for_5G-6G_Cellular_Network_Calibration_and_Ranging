import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles
from enviroment import ALUEnv

# ============================================================
#                       BASE TEST
# ============================================================


class BaseALUTest:
    """Sets up the environment, connects the DUT, and starts the components."""

    def __init__(self, dut):
        self.dut = dut

        # 1. Instantiate the Environment
        self.env = ALUEnv()

        # 2. Connect the DUT to the Environment blocks
        self.env.connect(self.dut)

    def run(self):
        # Start the hardware clock in the background
        cocotb.start_soon(Clock(self.dut.clk, 10, units="ns").start())

        # Fork the Agent and Scoreboard run phases into the background!
        cocotb.start_soon(self.env.agent.run())
        cocotb.start_soon(self.env.scoreboard.run())

    async def wait_completion(self, additional_cycles=10):
        # Wait for the pipeline to clear out before reporting
        await ClockCycles(self.dut.clk, additional_cycles)
        self.env.report()

# ============================================================
#               COCOTB TOP-LEVEL RUNNERS
# ============================================================


@cocotb.test()
async def run_reset_test(dut):
    test = BaseALUTest(dut)
    test.run()

    # Trigger the Generator block to create and send reset transactions
    await test.env.agent.generator.generate_reset(4)
    await test.wait_completion(5)


@cocotb.test()
async def run_add_test(dut):
    test = BaseALUTest(dut)
    test.run()

    # Trigger the Generator block to create and send ADD transactions
    await test.env.agent.generator.generate_add(21)
    await test.wait_completion(15)


@cocotb.test()
async def run_xor_test(dut):
    test = BaseALUTest(dut)
    test.run()

    # Trigger the Generator block to create and send XOR transactions
    await test.env.agent.generator.generate_xor(21)
    await test.wait_completion(15)


@cocotb.test()
async def run_and_test(dut):
    test = BaseALUTest(dut)
    test.run()

    # Trigger the Generator block to create and send AND transactions
    await test.env.agent.generator.generate_and(21)
    await test.wait_completion(15)


@cocotb.test()
async def run_or_test(dut):
    test = BaseALUTest(dut)
    test.run()

    # Trigger the Generator block to create and send OR transactions
    await test.env.agent.generator.generate_or(21)
    await test.wait_completion(15)


@cocotb.test()
async def run_random_test_with_reset(dut):
    test = BaseALUTest(dut)
    test.run()

    # Trigger the Generator block to create and send Random/Reset transactions
    await test.env.agent.generator.generate_random_with_reset(25)
    await test.wait_completion(15)
