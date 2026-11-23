"""Synthetic evaluator acceptance checks, no benchmark results."""
import importlib.util
from pathlib import Path
import sqlite3
import tempfile
import unittest
import pandas as pd

spec = importlib.util.spec_from_file_location('eval_calib', Path(__file__).with_name('eval_calib.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class CalibrationTest(unittest.TestCase):
    def test_counts_and_empty_threshold(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            library = root / 'lib.pqp'
            with sqlite3.connect(library) as conn:
                conn.executescript('''CREATE TABLE PRECURSOR(ID,DECOY);
                CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID,PEPTIDE_ID);
                CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID,PROTEIN_ID);
                CREATE TABLE PROTEIN(ID,PROTEIN_ACCESSION);
                INSERT INTO PRECURSOR VALUES(1,0),(2,0);
                INSERT INTO PRECURSOR_PEPTIDE_MAPPING VALUES(1,1),(2,2);
                INSERT INTO PEPTIDE_PROTEIN_MAPPING VALUES(1,1),(2,2);
                INSERT INTO PROTEIN VALUES(1,'human'),(2,'shuffle_p_target');''')
            for arm in ['never', 'always']:
                for seed in range(5):
                    directory = root / f'X_{arm}_s{seed}'
                    directory.mkdir()
                    pd.DataFrame(dict(precursor_id=[1,1,2,99], is_decoy=[0,0,0,1],
                                      score=[3,2,1,4], q_value=[.005,.005,.01,.001])).to_csv(
                                          directory / 'run.precursors.tsv', sep='\t', index=False)
            result = module.evaluate(str(library), str(root), 'X')
            self.assertEqual(len(result), 90)
            one = result[result.alpha == .01]
            self.assertTrue((one.N == 2).all())
            self.assertTrue((one.E == 1).all())
            self.assertTrue((one.combined == 1).all())
            self.assertTrue(result[result.alpha == .001].lower.isna().all())
            self.assertTrue((result[result.alpha == .005].N == 1).all())


if __name__ == '__main__':
    unittest.main()
