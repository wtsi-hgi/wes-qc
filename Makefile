.PHONY: test

test: unit-test integration-test

# test-ut-one-step:
# 	cd tests/unit_tests && pytest $(test)

test-it-one-step:
	cd tests/integration_tests && pytest test_integration.py::IntegrationTests::$(test)

integration-test:
	cd tests/integration_tests && pytest

integration-test-coverage:
	cd tests/integration_tests && pytest --cov=../..

unit-test:
	cd tests/unit_tests && pytest

unit-test-coverage:
	cd tests/unit_tests && pytest --cov=../.. 