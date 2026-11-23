"""Fabricated fixtures only. Run: python3 -m unittest discover -s <calib> -p 'test_asymmetry.py'."""
import csv
import gzip
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

from audit_asymmetry import audit


class AsymmetryTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.lib = self.root / 'synthetic.pqp'
        with sqlite3.connect(self.lib) as connection:
            connection.executescript('''
                CREATE TABLE PRECURSOR(ID INTEGER, TRAML_ID TEXT, DECOY INTEGER, CHARGE INTEGER);
                CREATE TABLE PEPTIDE(ID INTEGER, MODIFIED_SEQUENCE TEXT);
                CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID INTEGER, PEPTIDE_ID INTEGER);
            ''')
            for pid, name, side, charge, seq in [
                    (1, 'a2', 0, 2, 'AAA'), (2, 'a3', 0, 3, 'AAA'),
                    (3, 'DECOY_a2', 1, 2, 'AXA'), (4, 'DECOY_a3', 1, 3, 'AYA'),
                    (5, 'orphan', 0, 2, 'BBB')]:
                connection.execute('INSERT INTO PRECURSOR VALUES(?,?,?,?)', (pid, name, side, charge))
                connection.execute('INSERT INTO PEPTIDE VALUES(?,?)', (pid, seq))
                connection.execute('INSERT INTO PRECURSOR_PEPTIDE_MAPPING VALUES(?,?)', (pid, pid))
        self.features = self.root / 'synthetic.features.tsv'
        self.header = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 's_n', 's_co']
        self.rows = [
            [0, 1, 0, 0, 10, -1e30, 1, 1, 1],
            [1, 1, 0, 0, 20, 0, 1, 1, 0],
            [2, 2, 0, 0, 13, 0, 1, 1, 1],
            [3, 3, 1, 0, 10, 0, 1, 0, 0],
        ]
        self.write_features()
        self.reference = self.root / 'reference.tsv'
        self.reference.write_text('precursor_id\tis_decoy\tscore\n1\t0\t0\n3\t1\t0\n2\t0\t2\n4\t1\t0\n')

    def write_features(self):
        with self.features.open('w') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(self.header)
            writer.writerows(self.rows)

    def test_counts_pairs_zero_and_null_proxy(self):
        result = audit(self.lib, self.features, self.reference)
        precursors = {p['precursor_id']: p for p in result['precursors']}
        self.assertEqual(precursors['1']['candidate_windows'], 2)  # Includes NOT_KEPT.
        self.assertEqual(precursors['1']['windows_with_coeluting_sibling'], 1)
        self.assertEqual(precursors['4']['candidate_windows'], 0)
        self.assertIsNone(precursors['4']['max_coeluting_siblings'])
        self.assertEqual(result['unpaired_precursor_ids'], ['5'])
        pairs = {p['pair_key']: p for p in result['own_decoy_pairs']}
        self.assertEqual(pairs['a2']['delta_library_siblings'], 1)
        self.assertEqual(pairs['a2']['delta_candidate_windows'], 1)
        self.assertEqual(pairs['a2']['reference_region'], 'low_reference_both')
        self.assertEqual(pairs['a3']['reference_region'], 'outside_low_reference')
        self.assertEqual(audit(self.lib, self.features)['own_decoy_pairs'][0]['reference_region'], 'not_requested')

    def test_gzip_and_missing_reference(self):
        compressed = self.root / 'synthetic.features.tsv.gz'
        with gzip.open(compressed, 'wt') as handle:
            handle.write(self.features.read_text())
        self.reference.write_text('precursor_id\tis_decoy\tscore\n1\t0\t0\n')
        result = audit(self.lib, compressed, self.reference)
        self.assertTrue(all(p['reference_region'] == 'missing_reference' for p in result['own_decoy_pairs']))

    def test_invalid_inputs(self):
        for change in ('duplicate', 'side', 'nonfinite', 'sibling'):
            with self.subTest(change=change):
                original = [row[:] for row in self.rows]
                if change == 'duplicate':
                    self.rows.append(self.rows[0][:])
                elif change == 'side':
                    self.rows[0][2] = 1
                elif change == 'nonfinite':
                    self.rows[0][4] = 'nan'
                else:
                    self.rows[0][7] = 99
                self.write_features()
                with self.assertRaises(ValueError):
                    audit(self.lib, self.features)
                self.rows = original
        self.features.write_text('precursor_id\tis_decoy\tscore\n1\t0\t0\n')
        with self.assertRaises(ValueError):
            audit(self.lib, self.features)

    def test_library_only_null_candidates(self):
        result = audit(self.lib)
        self.assertEqual(result['candidate_population'], 'unavailable_library_only')
        for precursor in result['precursors']:
            for field in ('candidate_windows', 'candidate_siblings', 'windows_with_coeluting_sibling', 'max_coeluting_siblings'):
                self.assertIsNone(precursor[field])
        for pair in result['own_decoy_pairs']:
            self.assertEqual(pair['delta_library_siblings'], 1)
            self.assertIsNone(pair['delta_candidate_windows'])
        output = self.root / 'library-only.json'
        command = [sys.executable, str(Path(__file__).with_name('audit_asymmetry.py')),
                   '--library', str(self.lib), '--library-only', '--output', str(output)]
        self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
        subprocess.run(command + ['--synthetic'], check=True)
        self.assertIsNone(json.loads(output.read_text())['precursors'][0]['candidate_windows'])
        self.assertNotEqual(subprocess.run(command + ['--synthetic', '--reference', str(self.reference)], capture_output=True).returncode, 0)

    def test_cli_gate_and_synthetic_provenance(self):
        output = self.root / 'result.json'
        command = [sys.executable, str(Path(__file__).with_name('audit_asymmetry.py')),
                   '--library', str(self.lib), '--features', str(self.features), '--output', str(output)]
        blocked = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(blocked.returncode, 0)
        self.assertFalse(output.exists())
        subprocess.run(command + ['--synthetic'], check=True)
        result = json.loads(output.read_text())
        self.assertTrue(result['provenance']['synthetic'])
        self.assertEqual(len(result['provenance']['inputs']), 2)
        self.assertNotEqual(subprocess.run(command + ['--synthetic'], capture_output=True).returncode, 0)


if __name__ == '__main__':
    unittest.main()
