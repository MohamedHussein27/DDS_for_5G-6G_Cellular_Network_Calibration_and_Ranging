import cocotb
import logging

from cocotb.clock import Clock
from cocotb.triggers import RisingEdge
from cocotb.queue import Queue

from driver import ALUDriver
from monitor import ALUMonitor
from sequence import AddSequence, XorSequence, AndSequence, OrSequence
from scoreboard import ALUScoreboard
from coverage import ALUCoverage


# =====================================================
# Common Environment Setup
# =====================================================

async def base_test_setup(dut, sequence_class, count=500):

    cocotb.start_soon(Clock(dut.clk, 10, unit="ns").start())

    seq_queue = Queue()

    seq = sequence_class(seq_queue, count)

    drv = ALUDriver(dut, seq_queue)
    sb = ALUScoreboard()
    cov = ALUCoverage()
    mon = ALUMonitor(dut, sb, cov)

    await drv.reset()

    cocotb.start_soon(drv.run())
    cocotb.start_soon(mon.run())

    await seq.run()

    for _ in range(2000):
        await RisingEdge(dut.clk)

    logging.info("================================")
    logging.info(f"FINAL ERROR COUNT: {sb.error_count}")
    logging.info("================================")


# =====================================================
# ADD TEST
# =====================================================

@cocotb.test()
async def add_test(dut):

    logging.info("Starting ADD Test")

    await base_test_setup(dut, AddSequence)


# =====================================================
# XOR TEST
# =====================================================

@cocotb.test()
async def xor_test(dut):

    logging.info("Starting XOR Test")

    await base_test_setup(dut, XorSequence)


# =====================================================
# AND TEST
# =====================================================

@cocotb.test()
async def and_test(dut):

    logging.info("Starting AND Test")

    await base_test_setup(dut, AndSequence)


# =====================================================
# OR TEST
# =====================================================

@cocotb.test()
async def or_test(dut):

    logging.info("Starting OR Test")

    await base_test_setup(dut, OrSequence)