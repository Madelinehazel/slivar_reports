# Parse Slivar Reports

The scripts in this directory parse the tab-separated text file output from the [slivar](https://github.com/brentp/slivar) `tsv` function and add annotations to recreate the reports generate by [cre](https://github.com/ccmbioinfo/cre/blob/master/cre.sh)

## Requirements

Python3
- pandas
- argparse
- pytest

## Usage

```
python3 generate_report.py -report [slivar-generated tsv] -output [output csv]
```


## Tests

Tests for each function in parse_report/parse_functions.py are located in tests/test_parse_report.py. To install the parse_report package in editable mode, run:

```
pip install -e .
```
Now you can change the source code and run tests at will. To run tests:

```
python3 -m pytest
```

## Development

For each function your write, add tests to tests/test_parse_report.py. This allows us to thoroughly test code functionality for discrete features. 
