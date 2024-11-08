.PHONY: test

test:
	cd tests/unit_tests && pytest

coverage:
	cd tests/unit_tests && coverage run -m pytest && coverage report -m