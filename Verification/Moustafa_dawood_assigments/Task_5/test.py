import cocotb
import logging
from cocotb.clock import Clock
from cocotb.triggers import RisingEdge


from env import ALUEnv
from sequence import AddSequence, XorSequence, AndSequence, OrSequence

TOTAL_ERRORS = 0

async def base_test(dut, sequence):
    global TOTAL_ERRORS

    cocotb.start_soon(Clock(dut.clk,10,"ns").start())

    env = ALUEnv(dut)

    env.agent.seq = sequence

    env.connect()

    await env.run()

    for _ in range(1000):
        await RisingEdge(dut.clk)

    TOTAL_ERRORS += env.sb.error_count
   
    cocotb.log.info("========================================")
    cocotb.log.info(f"TEST COMPLETE: {sequence.__class__.__name__}")
    cocotb.log.info(f"Errors in THIS test: {env.sb.error_count}")
    cocotb.log.info(f" TOTAL ERRORS : {TOTAL_ERRORS}")
    cocotb.log.info("========================================")
    



@cocotb.test()
async def add_test(dut):

    seq = AddSequence(count=500)

    await base_test(dut, seq)


@cocotb.test()
async def xor_test(dut):

    seq = XorSequence(count=500)

    await base_test(dut, seq)


@cocotb.test()
async def and_test(dut):

    seq = AndSequence(count=500)

    await base_test(dut, seq)


@cocotb.test()
async def or_test(dut):

    seq = OrSequence(count=500)

    await base_test(dut, seq)