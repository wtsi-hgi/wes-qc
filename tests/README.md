# How to run tests

## Unit tests

Run all tests with: 
```Python
cd unit_tests/
python -m unittest test_all_modules
```

To run individual tests, you have to specify the test name. For instance:
```Python
python -m unittest test_all_modules.TestQCSteps.test_1_1_load_vcfs_to_mt
```

Unit tests can be executed in any order as they rely on reference test data stored in the S3 storage.

## Integration tests
Similarly to the unit tests, run all integration tests with:
```Python
cd integration_tests
python -m unittest test_integrations
```

Or select specific tests to run:
```Python
python -m unittest test_integrations.IntegrationTests.test_1_1_import_data
```
> [!NOTE]
> Integration tests test the whole pipeline steps for possible errors. Since some steps of the QC require the outputs of the previous steps, the execution order of the tests matters. Make sure to run the whole integration test suite first to get all the necessary files generated. Alternatively, you can opt for running separate tests, but keep in mind that this would require having all the inputs required by the tested step. 
