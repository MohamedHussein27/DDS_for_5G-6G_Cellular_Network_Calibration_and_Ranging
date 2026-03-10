import cocotb
from cocotb.clock import Clock
from cocotb.triggers import *
from environment import *

"""we will make three test cases"""

# 1. ADD XOR test
@cocotb.test()
async def add_xor_test(dut):
    

    # clock generation
    clk = Clock(dut.clk, 10, units="ns")
    await cocotb.start(clk.start())

    env = Environment()
    cocotb.log.info(" simulation starts ")

    env.connect()

    # variable to choose the sequence in the generator
    test_type = "add_xor_test"

    await env.run_environment(dut, test_type)

    await env.join_any.wait()


    # Report Phase
    cocotb.log.info("Completion signal received. Ending test.")
    env.s.report_test_cases() # s is the scoreboard instance

    # Print Coverage Results
    env.sub.report()

# 2. AND OR test
@cocotb.test()
async def and_or_test(dut):

    # clock generation
    clk = Clock(dut.clk, 10, units="ns")
    await cocotb.start(clk.start())

    # environment object
    env = Environment()
    cocotb.log.info(" simulation starts ")

    env.connect()

    # variable to choose the sequence in the generator
    test_type = "and_or_test"

    await env.run_environment(dut, test_type)

    await env.join_any.wait()


    # Report Phase
    cocotb.log.info("Completion signal received. Ending test.")
    env.s.report_test_cases() # s is the scoreboard instance

    # Print Coverage Results
    env.sub.report()

# 3. Random test
@cocotb.test()
async def random_test(dut):

    # clock generation
    clk = Clock(dut.clk, 10, units="ns")
    await cocotb.start(clk.start())

    # environment object
    env = Environment()
    cocotb.log.info(" simulation starts ")

    env.connect()

    # variable to choose the sequence in the generator
    test_type = "random_test"

    await env.run_environment(dut, test_type)

    await env.join_any.wait()


    # Report Phase
    cocotb.log.info("Completion signal received. Ending test.")
    env.s.report_test_cases() # s is the scoreboard instance

    # Print Coverage Results
    env.sub.report()
