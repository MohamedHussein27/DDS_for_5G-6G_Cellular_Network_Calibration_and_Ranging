import cocotb
from cocotb.clock import Clock
from cocotb.triggers import RisingEdge

from cocotb.queue import Queue

from driver import ALUDriver
from monitor import ALUMonitor
from sequence import ALUSequence
from scoreboard import ALUScoreboard
from coverage import ALUCoverage


@cocotb.test()
async def alu_test(dut):

    cocotb.start_soon(Clock(dut.clk, 10, unit="ns").start())

    seq_queue = Queue()

    seq = ALUSequence(seq_queue)
    drv = ALUDriver(dut, seq_queue)
    sb = ALUScoreboard()
    cov = ALUCoverage()
    mon = ALUMonitor(dut, sb, cov)

    await drv.reset()

    cocotb.start_soon(drv.run())
    cocotb.start_soon(mon.run())

    await seq.run()

    for _ in range(1050):
        await RisingEdge(dut.clk)

        # 3. Print the final results from  Scoreboard
    import logging
    logging.info("========================================")
    logging.info(f"FINAL ERROR COUNT: {sb.error_count}")
    logging.info("========================================")