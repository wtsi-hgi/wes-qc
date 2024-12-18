from utils.config import parse_config, parse_config_file, getp, parent_section, is_subsection, flatten, flat_to_nested, resolve_cvar_flat
from pathlib import PurePath
import unittest
import os

class TestConfigParser(unittest.TestCase):

    def test_parent_section(self):
        self.assertEqual(parent_section(''), '')
        self.assertEqual(parent_section('.'), '')
        self.assertEqual(parent_section('a'), '')
        self.assertEqual(parent_section('contra-spem-spero'), '')
        self.assertEqual(parent_section('a.b'), 'a')
        self.assertEqual(parent_section('contra-spem-spero.b'), 'contra-spem-spero')
        self.assertEqual(parent_section('a.b.c'), 'a.b')
        self.assertEqual(parent_section('contra.spem.spero'), 'contra.spem')

    def test_subsections(self):
        self.assertTrue(is_subsection('', ''))
        self.assertTrue(is_subsection('a', ''))
        self.assertTrue(is_subsection('contra-spem-spero', ''))
        self.assertTrue(is_subsection('a.b', 'a'))
        self.assertTrue(is_subsection('contra-spem-spero.b', 'contra-spem-spero'))
        self.assertTrue(is_subsection('a.b.c', 'a.b'))
        self.assertTrue(is_subsection('contra.spem.spero', 'contra.spem'))
        self.assertFalse(is_subsection('', 'a'))
        self.assertFalse(is_subsection('', 'contra-spem-spero'))
        self.assertFalse(is_subsection('a', 'a.b'))
        self.assertFalse(is_subsection('contra-spem-spero', 'contra-spem-spero.b'))
        self.assertFalse(is_subsection('a.b', 'a.b.c'))
        self.assertFalse(is_subsection('contra.spem', 'contra.spem.spero'))

    def test_flatten(self):
        d = {
            'a': {
                'b': 'a.b',
                'c': 'a.c',
                'd': { 'e': 'a.d.e', 'i1': 100 },
                'i2': 200,
                'f': { 'e': 'a.f.e', 'i4': 400 },
                'empty': dict()
            },
            'L1': [10, 20, 30],
            'L2': ['xx','yy','zz'],
            'b': dict(),
            'i3': 300
        }
        expected_d = {
            'a.b': 'a.b',
            'a.c': 'a.c',
            'a.d.e': 'a.d.e',
            'a.d.i1': 100,
            'a.i2': 200,
            'a.f.e': 'a.f.e',
            'a.f.i4': 400,
            # 'a.empty': dict(), # should not be added
            'L1': [10, 20, 30],
            'L2': ['xx','yy','zz'],
            # 'b': dict(), # should not be added
            'i3': 300
        }
        
        flat_d = flatten(d)
        
        self.assertDictEqual(expected_d, flat_d)
        
        # Order is important as well
        self.assertTupleEqual(tuple(expected_d.keys()), tuple(flat_d.keys()), "Order is wrong")

        # Change order
        expected_d.pop('a.b')
        expected_d['a.b'] = 'a.b'
        self.assertFalse(tuple(expected_d.keys()) == tuple(flat_d.keys()))

    def test_flat_to_nested(self):
        flat_d = {
            'a.b': 'a.b',
            'a.c': 'a.c',
            'a.d.e': 'a.d.e',
            'a.d.i1': 100,
            'a.i2': 200,
            'a.f.e': 'a.f.e',
            'a.f.i4': 400,
            # 'a.empty': dict(), # should not be added
            'L1': [10, 20, 30],
            'L2': ['xx','yy','zz'],
            # 'b': dict(), # should not be added
            'i3': 300
        }
        d = flat_to_nested(flat_d)
        flat_again_d = flatten(d, separator='$')
        expected_flat_again_d = { key.replace('.','$'): val for key,val in flat_d.items() }

        # Idempotence. Basically, flat_d == flatten(flat_to_nested(flat_d)).
        # Insertion/iteration order is important as well.
        self.assertTupleEqual(tuple(expected_flat_again_d.items()), tuple(flat_again_d.items()))

        # Ensure the dict is really nested
        self.assertEqual(d['a']['b'], 'a.b')
        self.assertEqual(d['a']['d']['e'], 'a.d.e')
        self.assertEqual(d['a']['d']['i1'], 100) # ensure types are ok

    def test_resolve_cvar(self):
        config = {
            'cvars.cvar0': 'a.x',
            'cvars.cvar1': 'x',
            'x': 'toplevel',
            'a.x': 'axxxx',
            'a.y': 'ayyyy',
            'b': '{cvar0}'
        }

        # config is FLAT on purpose
        self.assertEqual(resolve_cvar_flat('cvar0', config, 'field', True), 'axxxx') # cvar shortcut
        self.assertEqual(resolve_cvar_flat('a.x', config, 'field', True), 'axxxx') # full key
        self.assertEqual(resolve_cvar_flat('x', config, 'field', True), 'toplevel') # full key
        self.assertEqual(resolve_cvar_flat('x', config, 'a.other', True), 'axxxx') # local key in section a
        self.assertEqual(resolve_cvar_flat('.x', config, 'a.other', True), 'toplevel') # explicit full key / absolute key
        self.assertEqual(resolve_cvar_flat('cvar1', config, 'a.other', True), 'toplevel') # cvar shortcut points to the full key
        self.assertEqual(resolve_cvar_flat('a.nonexistant', config, '', True), '{a.nonexistant}')
        self.assertEqual(resolve_cvar_flat('a.nonexistant.nonexistant', config, '', True), '{a.nonexistant.nonexistant}')


    def test_cvars(self):
        expected_config = {
            'cvars.cvar_3': 'group_0.not_cvar',
            'cvars.cvar_custom': 'group_0.not_cvar', # This should be the second one, as custom cvars are inserted in the end, 
                                                     # and we are doing a depth-first iteration
            'cvar_1': '/dev/null/',
            'group_0.not_cvar': '/lorem/ipsum/',
            'group_1.cvar_2': '/foo/bar/',
            'group_2.f1_file': str(PurePath('/dev/null/baz.baz')),
            'group_2.f2_file': str(PurePath('/foo/bar/baz.baz')),
            'group_2.f3_file': str(PurePath('/lorem/ipsum/baz.baz')),
            'group_2.f4_file': str(PurePath('/lorem/ipsum/baz.baz')),
            'group_2.subgroup.f1_file': str(PurePath('/dev/null/baz.baz')),
            'group_2.subgroup.f2_file': str(PurePath('/foo/bar/baz.baz')),
            'group_2.subgroup.f3_file': str(PurePath('/lorem/ipsum/baz.baz')),
            'group_2.subgroup.f4_file': str(PurePath('/lorem/ipsum/baz.baz')),
            'group_2.a1': 100.0,
            'empty_file': '',
            'group_3.subgroup.f1': '12345',
            'group_3.subgroup.f2': '12345_1',
            'group_3.subgroup.f3': '12345_12',
            'group_3.f4': '12345_123',
            'group_3.f5': '12345_1234',
            'group_3.f1_1': '12345_1x12345_123',
        }
        config = parse_config('tests/config_tests/sample-config.yaml', { 'cvar_custom': 'group_0.not_cvar' })
        # print(config)        
        flat_config = flatten(config)
        self.assertTupleEqual(tuple(expected_config.items()), tuple(flat_config.items()))
        
        # ensure config is NESTED and not FLAT
        self.assertEqual(config['group_1']['cvar_2'], expected_config['group_1.cvar_2'])
        self.assertEqual(config['group_2']['f1_file'], expected_config['group_2.f1_file'])
        self.assertEqual(config['group_3']['subgroup']['f3'], expected_config['group_3.subgroup.f3'])

        # Test the whole file loading process

        os.environ['WES_CONFIG'] = os.path.abspath('tests/config_tests/sample-config.yaml')
        config = parse_config(additional_cvar_shortcuts={ 'cvar_custom': 'group_0.not_cvar' })
        flat_config = flatten(config)
        self.assertTupleEqual(tuple(expected_config.items()), tuple(flat_config.items()))

        # ensure config is NESTED and not FLAT
        self.assertEqual(config['group_1']['cvar_2'], expected_config['group_1.cvar_2'])
        self.assertEqual(config['group_2']['f1_file'], expected_config['group_2.f1_file'])
        self.assertEqual(config['group_3']['subgroup']['f3'], expected_config['group_3.subgroup.f3'])


    # TODO
    def test_undefined_cvar(self):
        with self.assertRaisesRegex(
            ValueError, 
            "cannot substitute undefined config variable 'b' referenced in field 'a'" 
            ) as cm:
            config = parse_config('tests/config_tests/undefined-cvar-config.yaml')

if __name__ == '__main__':
    unittest.main()

