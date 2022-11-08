# tfomics

A collection of functions and utilities for the tfomics project.

A python package is available in the gitlab package registry.
Please follow the instructions for [configuring pip to use the gitlab package registry](https://git.ecdf.ed.ac.uk/help/user/packages/pypi_repository/index.md) to install with pip.

To add the current dev version of this project to your poetry project, run

```bash
poetry add --git ssh://git@git.ecdf.ed.ac.uk/tfomics/prototypes/tfomics-utils.git
```

To install the current dev-version using pip, run
```bash
pip install --upgrade --force-reinstall git+ssh://git@git.ecdf.ed.ac.uk/tfomics/prototypes/tfomics-utils.git
```

# Code standards

Below is a wish list for code standards, at current no CI exists and so none of this is strongly enforced.

- Linting with pylint is enforced on all merges to master.
- We use pytest for unit tests and enforce 100% code coverage
- Code should be formatted with autopep8

## Unit tests

Unit tests can be found in [the test folder](https://git.ecdf.ed.ac.uk/tfomics/prototypes/tfomics-utils/tree/master/tests).

We use pytest for testing and pytest-cov for code coverage.
