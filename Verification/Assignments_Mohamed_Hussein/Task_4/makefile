PYTHON ?= python3
TEST = Testbench.py

run:
	@echo "Running test..."
	@$(PYTHON) $(TEST)

clean:
	find . -name "__pycache__" -type d -exec rm -r {} +
	find . -name "*.pyc" -delete
