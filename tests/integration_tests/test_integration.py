import inspect

from tests.integration_tests.integration_stub import IntegrationTestsStub

# /path/to/wes_qc must be in PYTHONPATH

PEDIGREE_FILE_PATH_TRIOS = """ '{metadir}/control_set_small_trios.ped' """


class IntegrationTests(IntegrationTestsStub):
    def __new__(cls, *args, **kwargs):
        # Find all stub methods in the parent class
        stub_methods = inspect.getmembers(
            IntegrationTestsStub, predicate=lambda x: inspect.isfunction(x) and x.__name__.startswith("stub_")
        )

        # For each stub method, create a corresponding test method
        for stub_name, stub_method in stub_methods:
            test_name = stub_name.replace("stub_", "test_")

            def create_test(stub):
                def test_method(self):
                    getattr(self, stub.__name__)()

                return test_method

            # Add the test method to the class
            setattr(cls, test_name, create_test(stub_method))

        return super().__new__(cls)

    @classmethod
    def setUpClass(cls):
        # Store pedigree file path as class variable
        cls.pedigree_file_path = PEDIGREE_FILE_PATH_TRIOS
        # Call parent's setUpClass without arguments
        super().setUpClass()
