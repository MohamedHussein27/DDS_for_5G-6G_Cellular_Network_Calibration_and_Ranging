import cocotb
import random
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles
from transaction import ALUTransaction
from enviroment import ALUEnv

# ============================================================
#               IMPORT YOUR SEQUENCE CLASSES
# ============================================================
from add_sequence import AddSequence
from xor_sequence import XorSequence
from and_sequence import AndSequence
from or_sequence import OrSequence
from random_sequence import RandomSequence

# ============================================================
#                       BASE TEST
# ============================================================


class BaseALUTest:
    """Sets up the environment, agent, scoreboard, and clock."""

    def __init__(self, dut):
        self.dut = dut
        self.env = ALUEnv(dut)

    def start(self):
        """Starts the clock and background components."""
        cocotb.start_soon(Clock(self.dut.clk, 10, units="ns").start())
        self.env.start()

    async def wait_completion(self):
        """Waits for transactions to pass through the DUT and reports."""
        await ClockCycles(self.dut.clk, 21)
        self.env.report()


# ============================================================
#                       TEST CLASSES
# ============================================================
class AddTest(BaseALUTest):
    async def run_test(self):
        self.start()
        seq = AddSequence()
        await seq.run(self.env.agent.sequencer)  # ROUTE TO SEQUENCER
        await self.wait_completion()


class XorTest(BaseALUTest):
    async def run_test(self):
        self.start()
        seq = XorSequence()
        await seq.run(self.env.agent.sequencer)  # ROUTE TO SEQUENCER
        await self.wait_completion()


class AndTest(BaseALUTest):
    async def run_test(self):
        self.start()
        seq = AndSequence()
        await seq.run(self.env.agent.sequencer)  # ROUTE TO SEQUENCER
        await self.wait_completion()


class OrTest(BaseALUTest):
    async def run_test(self):
        self.start()
        seq = OrSequence()
        await seq.run(self.env.agent.sequencer)  # ROUTE TO SEQUENCER
        await self.wait_completion()


class RandomTest(BaseALUTest):
    async def run_test(self):
        self.start()
        seq = RandomSequence()
        await seq.run(self.env.agent.sequencer)  # ROUTE TO SEQUENCER
        await self.wait_completion()

# ============================================================
#               COCOTB TOP-LEVEL RUNNERS
# ============================================================
# These small wrappers are required so cocotb can discover the tests.


@cocotb.test()
async def run_add_test(dut):
    test = AddTest(dut)
    await test.run_test()


@cocotb.test()
async def run_xor_test(dut):
    test = XorTest(dut)
    await test.run_test()


@cocotb.test()
async def run_and_test(dut):
    test = AndTest(dut)
    await test.run_test()


@cocotb.test()
async def run_or_test(dut):
    test = OrTest(dut)
    await test.run_test()


@cocotb.test()
async def run_random_test(dut):
    test = RandomTest(dut)
    await test.run_test()
