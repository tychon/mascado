try:
    from importlib.resources import read_text, contents, files  # type:ignore
except ModuleNotFoundError:
    from importlib_resources import read_text, contents, files  # type:ignore
from unittest import TestCase
import mascado.resources as res


class TestAccess (TestCase):
    def test_examples(self):
        affine_example = 'affine_example.txt'
        distorted_example = 'distorted_example.txt'
        resource_names = list(contents('mascado.resources.example_input'))

        assert affine_example in resource_names
        assert distorted_example in resource_names

        assert read_text('mascado.resources.example_input', affine_example)
        assert read_text('mascado.resources.example_input', distorted_example)
