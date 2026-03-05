import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles
from enviroment import ALUEnv

from add_sequence import AddSequence
from xor_sequence import XorSequence
from and_sequence import AndSequence
from or_sequence import OrSequence
from random_sequence import RandomSequence


class BaseALUTest:
    def __init__(self, dut):
        self.dut = dut
        self.env = ALUEnv(dut)

    def start_clock(self):
        cocotb.start_soon(Clock(self.dut.clk, 10, units="ns").start())

    async def setup_test(self):
        self.start_clock()
        self.env.start()
        await self.env.agent.driver.reset()

    async def wait_completion(self):
        await ClockCycles(self.dut.clk, 30)
        self.env.report()


class AddTest(BaseALUTest):
    async def run_test(self):
        await self.setup_test()
        seq = AddSequence()
        await seq.run(self.env.agent.sequencer)
        await self.wait_completion()


class XorTest(BaseALUTest):
    async def run_test(self):
        await self.setup_test()
        seq = XorSequence()
        await seq.run(self.env.agent.sequencer)
        await self.wait_completion()


class AndTest(BaseALUTest):
    async def run_test(self):
        await self.setup_test()
        seq = AndSequence()
        await seq.run(self.env.agent.sequencer)
        await self.wait_completion()


class OrTest(BaseALUTest):
    async def run_test(self):
        await self.setup_test()
        seq = OrSequence()
        await seq.run(self.env.agent.sequencer)
        await self.wait_completion()


class RandomTest(BaseALUTest):
    async def run_test(self):
        await self.setup_test()
        seq = RandomSequence()
        await seq.run(self.env.agent.sequencer)
        await self.wait_completion()


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
