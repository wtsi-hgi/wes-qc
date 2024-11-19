.PHONY: test

export PYTHONPATH:=$PYTHONPATH:$(shell pwd )
export PYSPARK_DRIVER_PYTHON:="/home/ubuntu/venv/bin/python"

test: unit-test integration-test

test-ut-one-step:
	cd tests/unit_tests && pytest -vv -s $(test)

test-it-one-step:
	cd tests/integration_tests && pytest -vv -s test_integration.py::IntegrationTests::$(test)

integration-test:
	cd tests/integration_tests && pytest

integration-test-coverage:
	cd tests/integration_tests && pytest --cov=../..

unit-test:
	cd tests/unit_tests && pytest

unit-test-coverage:
	cd tests/unit_tests && pytest --cov=../..

clear-logs:
	rm hail*.log
	rm hlrun_local_*
	rm tests/unit_tests/hail*.log
	rm tests/integration_tests/hail*.log
