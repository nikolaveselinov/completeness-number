PYTHON ?= python3
SAGE_PYTHON ?= ./run_sage.sh

.PHONY: help audit python-check shell-check native test-fast check test clean

help:
	@printf '%s\n' \
		'make audit      verify the published census without SageMath' \
		'make check      run all portable checks' \
		'make native     build the C++20 counters' \
		'make test-fast  run the Sage-free tests' \
		'make test       run the complete SageMath test suite' \
		'make clean      remove native build products'

audit:
	$(PYTHON) scripts/audit_census.py

python-check:
	$(PYTHON) -m compileall -q \
		certificates drinfeld_complete fast_prime fast_q2 scripts tests *.py

shell-check:
	bash -n run_sage.sh

native:
	$(MAKE) -C fast_q2 clean all
	$(MAKE) -C fast_prime clean all

test-fast:
	$(PYTHON) -m tests.test_release_census

check: audit python-check shell-check native test-fast

test:
	$(SAGE_PYTHON) -m pytest -q

clean:
	$(MAKE) -C fast_q2 clean
	$(MAKE) -C fast_prime clean
