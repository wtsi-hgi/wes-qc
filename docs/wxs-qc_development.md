# Developer's howto

This howto contains development howto and best practices.

## Set up a dev environment

Update your environment with extra dev dependencies:

```bash
uv sync
```

Set up pre-commit:

```bash
pre-commit install
```

This will set up the pre-commit hooks configured in `.pre-commit-config.yaml`, which include:
- Trailing whitespace removal
- End of file fixer
- YAML syntax checking
- Large file checking
- Ruff linting and formatting
- mypy type checking


## Run the tests and calculate coverage

The easiest way to run tests is using make utility and commands defined in `Makefile`.

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

## Use pre-commit hooks

Once `pre-commit` is installed and configured
(as described in the [Set up a dev environment](#set-up-a-dev-environment) section),
the hooks will automatically run on every commit.

### Run pre-commit manually

To run all pre-commit hooks on all files:
```bash
pre-commit run --all-files
```

To run pre-commit hooks on specific files:
```bash
pre-commit run --files <file1> <file2>
```

### Run MyPy type checking

`mypy` is configured to run in the manual stage because it currently produces many errors.
To run it:
```bash
pre-commit run mypy --hook-stage manual
```

### Skipping pre-commit hooks

In rare cases when you need to bypass pre-commit hooks (not recommended for regular use):
```bash
git commit -m "Your message" --no-verify
```
