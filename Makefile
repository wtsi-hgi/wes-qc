.PHONY: test

export PYTHONPATH:=$PYTHONPATH:$(shell pwd )
export PYSPARK_PYTHON:=/home/ubuntu/venv/bin/python
export PYSPARK_DRIVER_PYTHON:=/home/ubuntu/venv/bin/python

test: unit-test integration-test

test-ut-one-step:
	cd tests/unit_tests && pytest -vv -s -k $(test)

test-it-one-step:
	cd tests/integration_tests && pytest -vv -s -k $(test)

integration-test: clear-ht clear-logs
	cd tests/integration_tests && pytest


integration-test-coverage: clear-ht clear-logs
	cd tests/integration_tests && pytest --cov=../..


unit-test:
	cd tests/unit_tests && pytest

unit-test-coverage:
	cd tests/unit_tests && pytest --cov=../..

clear-logs:
	rm hail*.log || true
	rm hlrun_local_* || true
	rm tests/unit_tests/hail*.log || true
	rm tests/integration_tests/hail*.log || true

clear-ht:
	rm -rf tests/integration_tests/matrixtables/* || true


sync-to-private:
	git remote add origin-private git@github.com:wtsi-hgi/wes-qc-analysis.git || true
	git switch main
	git pull
	git push origin-private
