import cocotb
from cocotb.clock import Clock
from cocotb.triggers import Timer
from cocotb.queue import Queue  

from driver import Driver
from monitor import Monitor
from checker import Checker


@cocotb.test()
async def run_test(dut):

    # Clock
    clock = Clock(dut.clk, 10, units="ns")
    cocotb.start_soon(clock.start())

    # Reset
    dut.rst_n.value = 0
    dut.valid.value = 0
    await Timer(20, units="ns")
    dut.rst_n.value = 1

    # Queues
    in_queue = Queue()   
    out_queue = Queue()  

    # Components
    driver = Driver(dut, in_queue)
    monitor = Monitor(dut, out_queue)
    checker = Checker(in_queue, out_queue)

    # Start concurrent tasks
    cocotb.start_soon(driver.run())
    cocotb.start_soon(monitor.run())
    cocotb.start_soon(checker.run())

    await Timer(500, units="ns")
