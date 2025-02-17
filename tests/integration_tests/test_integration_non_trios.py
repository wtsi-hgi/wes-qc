import inspect

from tests.integration_tests.test_integration import IntegrationTestsStub

# TODO: this approach doesn't work. Leave as is. Find a netter way to parametrize tests: HSH-254
# INTEGRATION_TESTS_CONFIG_TEMPLATE = "config_test_template.yaml"
# INTEGRATION_TESTS_CONFIG_RENDERED_SAVEFILE = "integration_config_rendered_non_trios.yaml"

PEDIGREE_FILE_PATH_NON_TRIOS = "null"

# set up explicit runhash parameter for reproducible RF training
RF_RUN_TEST_HASH = "testhash_non-trios"  # manually set rf run id


class IntegrationTestsNonTrios(IntegrationTestsStub):
    def __new__(cls, *args, **kwargs):
        # Find all stub methods in the parent class
        stub_methods = inspect.getmembers(
            IntegrationTestsStub, predicate=lambda x: inspect.isfunction(x) and x.__name__.startswith("stub_")
        )

        # For each stub method, create a corresponding test method
        for stub_name, stub_method in stub_methods:
            test_name = stub_name.replace("stub_", "test_non_trios_")

            def create_test(stub):
                def test_method(self):
                    getattr(self, stub.__name__)()

                return test_method

            # Add the test method to the class
            setattr(cls, test_name, create_test(stub_method))

        return super().__new__(cls)

    @classmethod
    def setUpClass(cls):
        return super().setUpClass(pedifree_file_path=PEDIGREE_FILE_PATH_NON_TRIOS)
