.PHONY: test

test: unit-test integration-test

integration-test:
	cd tests/integration_tests && pytest

integration-test-coverage:
	cd tests/integration_tests && pytest --cov=../..

unit-test:
	cd tests/unit_tests && pytest

unit-test-coverage:
	cd tests/unit_tests && pytest --cov=../.. 