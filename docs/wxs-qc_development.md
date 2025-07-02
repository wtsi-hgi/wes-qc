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

## Development and code organization best practices

This section contains major suggestions sto maintain code style and structure.
Due to the limited number of developers,
the code refactoring is usually performed together with financial improvements.
Therefore, current guidelines represent the desired state for the codebase,
and not all pipeline parts follow it yet.

### Scripts organizations and sequence
- use numbered scripts for pipeline steps
- Break complex steps into smaller substeps using command-line arguments

### Main Function Structure
- Try to keep all data loading and saving (especially Hail structures) inside the main function.
- Copy existing fucntion structure for new scripts

### Pipeline step function design
- Accept matrix tables as primary input/output. Avoid writing/reading matrixtables inside the fucntion
- Use dictionary unpacking from parsed config for flexible argument passing
- Convert file paths to Spark format (`file://`) only before Hail operations

### Toolset Organization
- Maintain utility functions in the `wes-qc` folder
- Create and separate modules depending on the function purpose:
    - Data filterign/transformation
    - Statistical Analysis
    - Visualization

### Logging and reporting (TBD)
- Use Python's logging module to capture all important messages
- Provide informative error messages
- Consider making separate logs for
