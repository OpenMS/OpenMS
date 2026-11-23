"""Synthetic PQP tests, never benchmark libraries."""
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

from prepare_library import prepare, sha256


def fixture(path):
    with sqlite3.connect(path) as db:
        db.executescript('''CREATE TABLE PRECURSOR(ID INTEGER PRIMARY KEY, TRAML_ID TEXT, DECOY INTEGER, CHARGE INTEGER);
        CREATE TABLE PEPTIDE(ID INTEGER PRIMARY KEY, MODIFIED_SEQUENCE TEXT);
        CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID INTEGER, PEPTIDE_ID INTEGER);
        CREATE TABLE PROTEIN(ID INTEGER PRIMARY KEY, PROTEIN_ACCESSION TEXT);
        CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID INTEGER, PROTEIN_ID INTEGER);
        CREATE TABLE TRANSITION(ID INTEGER PRIMARY KEY, LIBRARY_INTENSITY REAL, PRODUCT_MZ REAL);
        CREATE TABLE TRANSITION_PRECURSOR_MAPPING(TRANSITION_ID INTEGER, PRECURSOR_ID INTEGER);
        CREATE TABLE TRANSITION_PEPTIDE_MAPPING(TRANSITION_ID INTEGER, PEPTIDE_ID INTEGER);
        INSERT INTO PROTEIN VALUES(1,'human'),(2,'ent_p_target');
        ''')
        for pid, name, decoy, sequence, size in [(1, 'a', 0, 'PEP', 20), (2, 'DECOY_a', 1, 'DPEP', 3),
                (3, 'b', 0, 'PEP', 0), (4, 'orphan', 0, 'ENT', 8), (5, 'DECOY_missing', 1, 'DENT', 1)]:
            db.execute('INSERT INTO PRECURSOR VALUES(?,?,?,2)', (pid, name, decoy))
            db.execute('INSERT INTO PEPTIDE VALUES(?,?)', (pid, sequence))
            db.execute('INSERT INTO PRECURSOR_PEPTIDE_MAPPING VALUES(?,?)', (pid, pid))
            db.execute('INSERT INTO PEPTIDE_PROTEIN_MAPPING VALUES(?,?)', (pid, 2 if pid == 4 else 1))
            for i in range(size):
                tid = pid * 100 + i
                db.execute('INSERT INTO TRANSITION VALUES(?,?,?)', (tid, i // 2, 300 + tid))
                db.execute('INSERT INTO TRANSITION_PRECURSOR_MAPPING VALUES(?,?)', (tid, pid))
                db.execute('INSERT INTO TRANSITION_PEPTIDE_MAPPING VALUES(?,?)', (tid, pid))


class PrepareTests(unittest.TestCase):
    def test_cli_preserves_source_metadata_short_and_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, search, scoring = [Path(tmp) / p for p in ('full.pqp', 'six.pqp', 'twelve.pqp')]
            fixture(source)
            before = sha256(source)
            result = subprocess.run([sys.executable, str(Path(__file__).with_name('prepare_library.py')),
                str(source), '--search6', str(search), '--scoring12', str(scoring)],
                check=True, capture_output=True, text=True, timeout=15)
            report = json.loads(result.stdout)
            self.assertEqual(before, sha256(source))
            self.assertTrue(report['source_unchanged'])
            self.assertEqual(report['source']['own_pairs'], 1)
            self.assertEqual(report['source']['own_pair_sibling_count_imbalance'], 1)
            self.assertEqual(report['source']['targets_without_own_decoy'], 2)
            self.assertEqual(report['source']['decoys_without_own_target'], 1)
            self.assertEqual(report['source']['entrapment_structural_metadata']['ratio'], .5)
            self.assertGreater(report['source']['numeric_id_intensity_order_violations'], 0)
            for output, limit in [(search, 6), (scoring, 12)]:
                with sqlite3.connect(output) as db, sqlite3.connect(source) as src:
                    for table in ('PRECURSOR', 'PEPTIDE', 'PROTEIN', 'PRECURSOR_PEPTIDE_MAPPING', 'PEPTIDE_PROTEIN_MAPPING'):
                        self.assertEqual(db.execute(f'SELECT * FROM {table}').fetchall(), src.execute(f'SELECT * FROM {table}').fetchall())
                    actual = db.execute('SELECT t.LIBRARY_INTENSITY,t.PRODUCT_MZ FROM TRANSITION t JOIN TRANSITION_PRECURSOR_MAPPING m ON m.TRANSITION_ID=t.ID WHERE m.PRECURSOR_ID=1 ORDER BY t.ID').fetchall()
                    expected = src.execute('SELECT LIBRARY_INTENSITY,PRODUCT_MZ FROM TRANSITION WHERE ID BETWEEN 100 AND 119 ORDER BY LIBRARY_INTENSITY DESC, ID LIMIT ?', (limit,)).fetchall()
                    self.assertEqual(actual, expected)
                    self.assertEqual(db.execute('SELECT count(*) FROM TRANSITION_PEPTIDE_MAPPING').fetchone(), db.execute('SELECT count(*) FROM TRANSITION').fetchone())
                    self.assertEqual(db.execute('SELECT count(*) FROM TRANSITION_PRECURSOR_MAPPING WHERE PRECURSOR_ID=2').fetchone()[0], 3)
            with sqlite3.connect(search) as six, sqlite3.connect(scoring) as twelve:
                scoring_rows = dict((r[0], r) for r in twelve.execute('SELECT * FROM TRANSITION'))
                for row in six.execute('SELECT * FROM TRANSITION'):
                    self.assertEqual(row, scoring_rows[row[0]])
            with self.assertRaises(ValueError):
                prepare(source, search, scoring)
            self.assertEqual(before, sha256(source))

    def test_rejects_mixed_labels_and_inconsistent_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, search, scoring = [Path(tmp) / p for p in ('full.pqp', 'six.pqp', 'twelve.pqp')]
            fixture(source)
            with sqlite3.connect(source) as db:
                target = db.execute('SELECT ID FROM PRECURSOR WHERE DECOY=0 LIMIT 1').fetchone()[0]
                peptide = db.execute('SELECT PEPTIDE_ID FROM PRECURSOR_PEPTIDE_MAPPING WHERE PRECURSOR_ID=?', (target,)).fetchone()[0]
                db.execute("INSERT INTO PROTEIN VALUES(999,'extra_p_target')")
                db.execute('INSERT INTO PEPTIDE_PROTEIN_MAPPING VALUES(?,999)', (peptide,))
            with self.assertRaisesRegex(ValueError, 'mixed target'):
                prepare(source, search, scoring)
            with sqlite3.connect(source) as db:
                db.execute('DELETE FROM PEPTIDE_PROTEIN_MAPPING WHERE PROTEIN_ID=999')
                db.execute('ALTER TABLE PEPTIDE ADD COLUMN DECOY INTEGER DEFAULT 0')
            with self.assertRaisesRegex(ValueError, 'peptide decoy flags'):
                prepare(source, search, scoring)
            self.assertFalse(search.exists())

    def test_shared_target_decoy_peptide_from_decoy_generator(self):
        """LABEL_AMENDMENT: a decoy equal to a target sequence shares its PEPTIDE row and adds a DECOY_ protein."""
        with tempfile.TemporaryDirectory() as tmp:
            source, search, scoring = [Path(tmp) / p for p in ('full.pqp', 'six.pqp', 'twelve.pqp')]
            fixture(source)
            with sqlite3.connect(source) as db:
                db.execute('ALTER TABLE PEPTIDE ADD COLUMN DECOY INTEGER DEFAULT 0')
                db.execute('UPDATE PEPTIDE SET DECOY=1 WHERE ID IN (2,5)')
                # Entrapment target 4 (ENT) and decoy 2 now share peptide 4, flagged as decoy, with a DECOY_ protein.
                db.execute("INSERT INTO PROTEIN VALUES(998,'DECOY_human')")
                db.execute('UPDATE PRECURSOR_PEPTIDE_MAPPING SET PEPTIDE_ID=4 WHERE PRECURSOR_ID=2')
                db.execute('UPDATE PEPTIDE SET DECOY=1 WHERE ID=4')
                db.execute('INSERT INTO PEPTIDE_PROTEIN_MAPPING VALUES(4,998)')
            report = prepare(source, search, scoring)
            shared = report['source']['target_decoy_shared_peptides']
            self.assertEqual((shared['peptides'], shared['target_precursors'], shared['decoy_precursors']), (1, 1, 1))
            # The DECOY_ protein does not change the entrapment label of target 4 or the ratio.
            self.assertEqual(report['source']['entrapment_structural_metadata']['ratio'], .5)
            with sqlite3.connect(source) as db:
                db.execute('UPDATE PEPTIDE SET DECOY=1 WHERE ID=1')
            for p in (search, scoring):
                p.unlink()
            with self.assertRaisesRegex(ValueError, 'peptide decoy flags'):
                prepare(source, search, scoring)
            with sqlite3.connect(source) as db:
                db.execute('UPDATE PEPTIDE SET DECOY=0 WHERE ID=1')
                db.execute('DELETE FROM PEPTIDE_PROTEIN_MAPPING WHERE PEPTIDE_ID=4 AND PROTEIN_ID=2')
            with self.assertRaisesRegex(ValueError, 'mixed target'):
                prepare(source, search, scoring)

    def test_rejects_excess_fragments_and_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, search, scoring = [Path(tmp) / p for p in ('full.pqp', 'six.pqp', 'twelve.pqp')]
            fixture(source)
            with self.assertRaises(ValueError):
                prepare(source, source, scoring)
            with sqlite3.connect(source) as db:
                db.execute('INSERT INTO TRANSITION VALUES(999,1,500)')
                db.execute('INSERT INTO TRANSITION_PRECURSOR_MAPPING VALUES(999,1)')
            with self.assertRaisesRegex(ValueError, 'exceeds 20'):
                prepare(source, search, scoring)
            self.assertFalse(search.exists())


if __name__ == '__main__':
    unittest.main()
