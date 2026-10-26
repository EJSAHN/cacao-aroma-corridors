"""Input/output protection and stage orchestration tests."""
import argparse
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('revision_runner', ROOT / 'run_revision.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class InterfaceTests(unittest.TestCase):
    def context(self, root):
        return runner.Context('background', root / 'supplement.xlsx', root / 'sources', root / 'output')

    def test_altered_source_fails_hash_check(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path = root / 'source'
            path.write_bytes(b'original')
            expected = runner.digest(path)
            context = self.context(root)
            context.verify_file(path, expected, 'test')
            path.write_bytes(b'altered')
            with self.assertRaises(ValueError):
                context.verify_file(path, expected, 'test')

    def test_output_cannot_escape_table_directory(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(ValueError):
                self.context(Path(temp)).write_csv('../../escape.csv', [{'count': 1}])

    def test_source_output_overlap_rejected_before_mutation(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args = argparse.Namespace(stage='background', supplement_s1=str(root / 'supplement.xlsx'),
                                      source_root=str(root / 'sources'), output_root=str(root))
            with self.assertRaises(ValueError):
                runner.run_one(args)
            self.assertEqual(list(root.iterdir()), [])

    def test_csv_null_and_structured_values_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.context(root).write_csv('values.csv', [{'missing': None, 'ids': {'b', 'a'}, 'zero': 0}])
            with (root / 'output/tables/values.csv').open(encoding='utf-8-sig') as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual(row['missing'], '')
            self.assertEqual(json.loads(row['ids']), ['a', 'b'])
            self.assertEqual(row['zero'], '0')

    def test_all_stages_stop_after_failure(self):
        args = ['run_revision.py', '--stage', 'all', '--supplement-s1', 'S1.xlsx',
                '--source-root', 'sources', '--output-root', 'results']
        with patch('sys.argv', args), patch.object(runner.subprocess, 'run') as execute:
            execute.side_effect = [argparse.Namespace(returncode=0), argparse.Namespace(returncode=2)]
            self.assertEqual(runner.main(), 2)
            self.assertEqual(execute.call_count, 2)
            stages = [call.args[0][call.args[0].index('--stage') + 1] for call in execute.call_args_list]
            self.assertEqual(stages, ['harmonization', 'genotypes'])

    def test_scientific_modules_match_recorded_bytes(self):
        records = json.loads((ROOT / 'SCIENTIFIC_SOURCE_MANIFEST.json').read_text())
        self.assertTrue(records)
        for row in records:
            self.assertEqual(runner.digest(ROOT / row['file']), row['sha256'], row['file'])


if __name__ == '__main__':
    unittest.main()
