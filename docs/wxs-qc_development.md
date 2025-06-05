# Developer's howto

## Develop dependencies

For developing the pipeline code, install extra dev and test dependencies:

```bash
uv pip install -r pyproject.toml --extra dev
uv pip install -r pyproject.toml --extra test
```


## How to run the tests and calculate coverage

The tests currently require running on the SPARK cluster. There are plans to make them runnable locally.
They can be run by commands defined in `Makefile`.

To run all the tests:
```bash
make test
```
Or you can specify the type of test to run
```bash
make unit-test
make integration-test
```

To run the tests with coverage:
```bash
make unit-test-coverage
make integration-test-coverage
```

## To run pre-commit hooks on commit

1. Install pre-commit
```shell
pip install pre-commit
```
2. `pre-commit` will automatically run on every commit
3. To run pre-commit manually on specific files
```shell
pre-commit run --files <file1> <file2>
```
4. `mypy` is configured to run manually because now it produces too many errors. To run it:
```shell
pre-commit run --hook-stage manual
```
