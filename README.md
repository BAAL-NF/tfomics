# tfomics

A collection of functions and utilities for the tfomics project.

To add this project to your poetry project, run

```bash
poetry add --git ssh://git@git.ecdf.ed.ac.uk/oalmelid/tfomics.git
```

# Code standards

Below is a wish list for code standards, at current no CI exists and so none of this is strongly enforced.

- Linting with pylint is enforced on all merges to master.
- We use pytest for unit tests and enforce 100% code coverage
- Code should be formatted with autopep8

## Unit tests

Unit tests can be found in [the test folder](https://git.ecdf.ed.ac.uk/oalmelid/tfomics/tree/master/tests).

We use pytest for testing and pytest-cov for code coverage.
