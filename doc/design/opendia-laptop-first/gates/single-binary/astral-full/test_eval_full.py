import importlib.util
from pathlib import Path
import sqlite3
import tempfile
import unittest
import subprocess
import sys
import pandas as pd

spec=importlib.util.spec_from_file_location('eval_full',Path(__file__).with_name('eval_full.py'))
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

class FullEvaluationTest(unittest.TestCase):
    def test_failures_counts_and_resources(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            lib=root/'lib.pqp'
            with sqlite3.connect(lib) as conn:
                conn.executescript('''CREATE TABLE PRECURSOR(ID,DECOY);
                CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID,PEPTIDE_ID);
                CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID,PROTEIN_ID);
                CREATE TABLE PROTEIN(ID,PROTEIN_ACCESSION);
                INSERT INTO PRECURSOR VALUES(1,0),(2,0);
                INSERT INTO PRECURSOR_PEPTIDE_MAPPING VALUES(1,1),(2,2);
                INSERT INTO PEPTIDE_PROTEIN_MAPPING VALUES(1,1),(2,2);
                INSERT INTO PROTEIN VALUES(1,'human'),(2,'shuffle_p_target');''')
            run=root/'raw01_s0'
            run.mkdir()
            (run/'exit_code.txt').write_text('0\n')
            (run/'time.txt').write_text('Maximum resident set size (kbytes): 1024\nElapsed (wall clock) time (h:mm:ss or m:ss): 1:02.00\n')
            pd.DataFrame(dict(precursor_id=[1,2,3],is_decoy=[0,0,1],score=[3,2,4],
                              q_value=[.005,.01,.001],sequence=['ABC','DEF','GHI'])).to_csv(
                run/'LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv',sep='\t',index=False)
            counts,status=module.evaluate(lib,root)
            self.assertEqual(len(status),15)
            self.assertEqual((status.status=='SUCCESS').sum(),1)
            self.assertEqual(len(counts),9)
            row=counts[counts.alpha==.01].iloc[0]
            self.assertEqual((row.N,row.E,row.human_precursors),(2,1,1))
            self.assertEqual(row.combined,1)
            self.assertTrue(pd.isna(counts[counts.alpha==.001].iloc[0].lower))
            self.assertTrue(status.iloc[0].within_32GiB_observed_rss)
            self.assertEqual(status.iloc[0].wall_time,'1:02.00')
            file = run/'LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv'
            frame = pd.read_csv(file, sep='\t')
            pd.concat([frame, frame.iloc[[0]]]).to_csv(file, sep='\t', index=False)
            with self.assertRaisesRegex(ValueError, 'duplicate target winners'):
                module.evaluate(lib, root)

    def test_literal_case_sensitive_entrapment_labels(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            lib = root/'lib.pqp'
            # Include LIKE wildcard/case false hits, a literal substring, mixed
            # mappings, and decoy mappings that must not affect target labels.
            accessions = [
                'sp|Q123P_target',
                'shuffle_P_TARGET',
                'shuffle_p_Target',
                'shuffle_pXtarget',
                'prefix_p_target_suffix',
                'human',
                'DECOY_shuffle_p_target',
                'DECOY_human',
            ]
            with sqlite3.connect(lib) as conn:
                conn.executescript('''CREATE TABLE PRECURSOR(ID,DECOY);
                CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID,PEPTIDE_ID);
                CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID,PROTEIN_ID);
                CREATE TABLE PROTEIN(ID,PROTEIN_ACCESSION);''')
                conn.executemany('INSERT INTO PROTEIN VALUES(?,?)', enumerate(accessions, 1))
                conn.executemany('INSERT INTO PRECURSOR VALUES(?,0)', [(i,) for i in range(1, 8)])
                conn.executemany('INSERT INTO PRECURSOR_PEPTIDE_MAPPING VALUES(?,?)',
                                 [(i, i) for i in range(1, 8)])
                conn.executemany('INSERT INTO PEPTIDE_PROTEIN_MAPPING VALUES(?,?)',
                                 [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5),
                                  (6, 6), (6, 7), (7, 5), (7, 6), (7, 8)])
            run = root/'raw01_s0'
            run.mkdir()
            (run/'exit_code.txt').write_text('0\n')
            # Distinct cutoffs expose each individual label, not just totals.
            pd.DataFrame(dict(precursor_id=list(range(1, 8)), is_decoy=[0]*7,
                              score=[3]*7, q_value=module.ALPHAS[:7],
                              sequence=[f'PEPTIDE{i}' for i in range(1, 8)])).to_csv(
                run/'LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv', sep='\t', index=False)
            counts, _ = module.evaluate(lib, root)
            self.assertEqual(counts.E.tolist(), [0, 0, 0, 0, 1, 1, 2, 2, 2])
            self.assertEqual(counts.human_precursors.tolist(), [1, 2, 3, 4, 4, 5, 5, 5, 5])
            self.assertTrue((counts.r == 2/5).all())

    def test_curve_plot_synthetic(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = [dict(raw=raw, alpha=q, lower_median=q/2, lower_min=q/3,
                         lower_max=q*.8, combined_median=q, combined_min=q*.8,
                         combined_max=q*1.2)
                    for raw in ['01','02','03'] for q in [.001,.01,.05]]
            pd.DataFrame(rows).to_csv(root/'calibration_curve.tsv', sep='\t', index=False)
            subprocess.run([sys.executable, str(Path(__file__).with_name('plot_full.py')), str(root)], check=True)
            self.assertGreater((root/'calibration_curves.png').stat().st_size, 100)
            self.assertGreater((root/'calibration_curves.svg').stat().st_size, 100)

if __name__=='__main__':
    unittest.main()
