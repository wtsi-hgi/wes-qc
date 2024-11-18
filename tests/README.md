# How to run tests

## Regression tests

Run all tests with: 
```Python
cd regression_tests/
python -m unittest test_regression
```

To run individual tests, you have to specify the test name. For instance:
```Python
python -m unittest test_regression.RegressionTests.test_1_1_load_vcfs_to_mt
```

Regression tests can be executed in any order as they rely on the reference test data stored in a dedicated S3 bucket.

## Integration tests
Before running integration tests, make sure that the `PYSPARK_PYTHON` environment variable is set correctly. 

Similarly to regression tests, run all integration tests with:
```Python
cd integration_tests
python -m unittest test_integration
```

Or select specific tests to run:
```Python
python -m unittest test_integration.IntegrationTests.test_1_1_import_data
```
> [!NOTE]
> Integration tests check all the pipeline steps for possible errors. Since some steps of the QC require the outputs of the previous ones, the execution order of the tests matters. Make sure to run the whole integration test suite first to get all the necessary files generated if you are going to run individual integration tests later.
