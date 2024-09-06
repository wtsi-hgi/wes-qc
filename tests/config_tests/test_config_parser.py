from utils.config import parse_config, getp, parent_section, is_subsection, flatten, resolve_cvar
from pathlib import PurePath
import unittest

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


    def test_resolve_cvar(self):
        config = {
            'cvar_shortcuts.cvar0': 'a.x',
            'cvar_shortcuts.cvar1': 'x',
            'x': 'toplevel',
            'a.x': 'axxxx',
            'a.y': 'ayyyy',
            'b': '{cvar0}'
        }
        self.assertEqual(resolve_cvar('cvar0', config, 'field', True), 'axxxx') # cvar shortcut
        self.assertEqual(resolve_cvar('a.x', config, 'field', True), 'axxxx') # full key
        self.assertEqual(resolve_cvar('x', config, 'field', True), 'toplevel') # full key
        self.assertEqual(resolve_cvar('x', config, 'a.other', True), 'axxxx') # local key in section a
        self.assertEqual(resolve_cvar('.x', config, 'a.other', True), 'toplevel') # explicit full key / absolute key
        self.assertEqual(resolve_cvar('cvar1', config, 'a.other', True), 'toplevel') # cvar shortcut points to the full key


    def test_cvars(self):
        expected_config = {
            'cvar_shortcuts.cvar_custom': 'group_0.not_cvar',
            'cvar_shortcuts.cvar_3': 'group_0.not_cvar',
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
        print(config)        

        for keyp, correct_val in expected_config.items():
            val = config[keyp]
            self.assertEqual(val, correct_val, f"config[{keyp}] is '{val}' ({type(val)}) and not '{correct_val}' ({type(val)})")
        self.assertEqual(len(config), len(expected_config), str(set(config.keys()).symmetric_difference(expected_config.keys())))

    # TODO
    def test_undefined_cvar(self):
        with self.assertRaisesRegex(
            ValueError, 
            "cannot substitute undefined config variable 'b' referenced in field 'a'" 
            ) as cm:
            config = parse_config('tests/config_tests/undefined-cvar-config.yaml')


if __name__ == '__main__':
    unittest.main()

