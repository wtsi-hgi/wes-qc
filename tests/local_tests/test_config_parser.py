from utils.config import parse_config, getp
from pathlib import PurePath
import unittest

# _file to activate path mode
correct_parsed_config = {
    'cvar_1': '/dev/null/',
    'group_0.not_cvar': '/lorem/ipsum/',
    'group_1.cvar_2': '/foo/bar/',
    'group_2.f1_file': str(PurePath('/dev/null/baz.baz')),
    'group_2.f2_file': str(PurePath('/foo/bar/baz.baz')),
    'group_2.f3_file': str(PurePath('/lorem/ipsum/baz.baz')),
    'group_2.subgroup.f1_file': str(PurePath('/dev/null/baz.baz')),
    'group_2.subgroup.f2_file': str(PurePath('/foo/bar/baz.baz')),
    'group_2.subgroup.f3_file': str(PurePath('/lorem/ipsum/baz.baz')),
    'group_2.a1': 100.0
}

class TestConfigParser(unittest.TestCase):

    def test_cvars(self):
        config = parse_config('tests/local_tests/sample-config.yaml', { 'cvar_3': 'group_0.not_cvar' })
        
        for keyp, correct_val in correct_parsed_config.items():
            val = getp(config, keyp)
            self.assertEqual(val, correct_val, f"config[{keyp}] is '{val}' ({type(val)}) and not '{correct_val}' ({type(val)})")
    


if __name__ == '__main__':
    unittest.main()
