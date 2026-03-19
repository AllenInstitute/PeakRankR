# Contributing to PyPeakRanker

Thank you for your interest in contributing!

## Reporting issues

Please use the [GitHub issue tracker](https://github.com/AllenInstitute/PyPeakRankR/issues).
Include:
- A minimal reproducible example
- Your Python version and OS
- The full error message and traceback

## Development setup

```bash
git clone https://github.com/AllenInstitute/PyPeakRankR
cd PyPeakRankR
pip install -e ".[dev]"
```

## Running tests

```bash
pytest tests/
```

## Code style

This project uses `black` for formatting and `ruff` for linting.

```bash
black src/ tests/
ruff check src/ tests/
```

## Pull requests

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Commit your changes with clear messages
4. Open a pull request against `main`

Please ensure all tests pass and new features include tests.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
