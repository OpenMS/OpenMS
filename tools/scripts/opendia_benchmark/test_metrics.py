#!/usr/bin/env python3
"""Regression tests for strict benchmark views, zero runs and native quantification."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from benchmark_metrics import deduplicate, entity_metrics, mapping_kind, stats, views
from validate_timsim import best_results, read_tsv


class MetricsTest(unittest.TestCase):
    def row(self, **changes):
        row = dict(run_name='r1', Sequence='PEPTIDE', FullPeptideName='PEPTIDE', Charge='2',
                   ProteinName='P', Intensity='10', m_score='.001',
                   m_score_peptide_global='.001', m_score_protein_global='.001')
        row.update(changes)
        return row

    def view(self, rows, path=Path('/missing/OpenDIA.results.tsv')):
        return views(rows, ['r1', 'r2', 'r3'], {r: 'replicate' for r in ['r1', 'r2', 'r3']}, path, read_tsv)

    def test_duplicate_low_q_then_high_intensity(self):
        rows = [self.row(Intensity='100', m_score='.005'), self.row(), self.row(Intensity='20')]
        self.assertEqual(next(iter(deduplicate(rows).values()))['Intensity'], '20')
        self.assertEqual(len(deduplicate(rows + [self.row(FullPeptideName='PEP(UniMod:1)TIDE')])), 2)

    def test_no_precursor_fallback(self):
        for invalid in ['', 'nan', 'inf', '-.1', '.02']:
            self.assertEqual(deduplicate([self.row(m_score=invalid, ms2_m_score='0')]), {})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'results.tsv'
            path.write_text('run_name\tSequence\tCharge\tms2_m_score\nr1\tPEPTIDE\t2\t0\n')
            with self.assertRaises(SystemExit):
                best_results(path, .01, ('ms2_m_score',))

    def test_global_filters_are_additional(self):
        result = self.view([self.row(m_score='.02'), self.row(Sequence='OTHER', FullPeptideName='OTHER', m_score_peptide_global='.02')])
        self.assertEqual(result['precursor_only']['per_run']['r1']['precursor'], 1)
        self.assertEqual(result['standard']['status'], 'ok')
        self.assertEqual(result['standard']['per_run']['r1']['precursor'], 0)
        self.assertEqual(result['standard']['per_run']['r3']['precursor'], 0)
        self.assertEqual(result['precursor_only']['count_summary']['replicate']['precursor']['values'], [1, 0, 0])

    def test_missing_nonfinite_global_is_na(self):
        for value in ['', 'nan', 'inf']:
            self.assertEqual(self.view([self.row(m_score_protein_global=value)])['standard']['status'], 'NA')
        row = self.row()
        del row['m_score_peptide_global']
        self.assertEqual(self.view([row])['standard']['status'], 'NA')

    def test_discarded_rows_do_not_disable_global_view(self):
        losing = self.row(m_score='.005', m_score_protein_global='nan')
        outside = self.row(FullPeptideName='OTHER', m_score='.04', m_score_peptide_global='')
        result = self.view([self.row(), losing, outside])
        self.assertEqual(result['standard']['status'], 'ok')
        self.assertEqual(result['standard']['per_run']['r1']['precursor'], 1)
        empty = self.view([outside])
        self.assertEqual(empty['standard']['status'], 'ok')
        self.assertEqual(empty['standard']['per_run']['r1']['precursor'], 0)
        missing_schema = dict(outside)
        del missing_schema['m_score_peptide_global']
        self.assertEqual(self.view([missing_schema])['standard']['status'], 'NA')

    def test_sample_sd_and_correlation_support(self):
        self.assertEqual(stats([1, 2, 3])['sample_sd'], 1)
        metrics = entity_metrics({'a': {'r1': 1, 'r2': 2, 'r3': 3}, 'b': {'r1': 4}}, ['r1', 'r2', 'r3'])
        self.assertAlmostEqual(metrics['cv_values'][0], .5)
        self.assertEqual(metrics['cv_n'], 1)
        self.assertEqual(metrics['pairwise_log2_correlations']['r1__r2']['n'], 1)
        self.assertIsNone(metrics['pairwise_log2_correlations']['r1__r2']['pearson_log2'])
        self.assertEqual(entity_metrics({}, ['r1', 'r2'])['complete_entities'], 0)

    def test_native_matrix_not_python_sum(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'OpenDIA.results.tsv'
            for level, header in [('precursor', 'FullPeptideName\tCharge\tProteinName'), ('peptide', 'FullPeptideName'), ('protein', 'ProteinName')]:
                identity = 'PEPTIDE\t2\tP' if level == 'precursor' else ('P' if level == 'protein' else 'PEPTIDE')
                path.with_name(f'OpenDIA.{level}.matrix.tsv').write_text(f'{header}\tr1\tr2\tr3\n{identity}\t2\t4\t6\n')
            view = self.view([self.row(run_name=r, Intensity='999') for r in ['r1', 'r2', 'r3']], path)
            quant = view['standard']['quantification']
            self.assertEqual(quant['protein_label']['source'], 'native_matrix')
            self.assertAlmostEqual(quant['protein_label']['conditions']['replicate']['cv_values'][0], .5)
            failed = self.view([self.row(m_score_peptide_global='.02')], path)
            self.assertEqual(failed['standard']['quantification']['precursor']['conditions']['replicate']['entities'], 0)

    def test_astral_cli_keeps_zero_runs_and_excludes_mixed_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'results.tsv'
            rows = [self.row(), self.row(run_name='r2', m_score='.02'),
                    self.row(FullPeptideName='OTHER', ProteinName='A_target;B_p_target')]
            import csv
            with path.open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
                writer.writeheader()
                writer.writerows(rows)
            completed = subprocess.run([sys.executable, str(Path(__file__).with_name('metrics_astral.py')), str(path), '--runs', 'r1', 'r2', 'r3'], check=True, capture_output=True, text=True)
            result = json.loads(completed.stdout)
            self.assertEqual(result['runs'], ['r1', 'r2', 'r3'])
            self.assertEqual(result['views']['standard']['per_run']['r1']['precursor'], 1)
            self.assertEqual(result['views']['standard']['per_run']['r3']['precursor'], 0)
            self.assertEqual(result['identification']['ambiguous_mixed_mapping_rows'], 1)
            self.assertEqual(result['fdr_curve_per_run']['r3'][-1]['rows'], 0)

    def test_fixed_cli_thresholds(self):
        for script, paths in [('metrics_astral.py', ['missing.tsv']), ('validate_timsim.py', ['dataset', 'target.tsv', 'entrapment.tsv'])]:
            for value in ['.005', '.1']:
                completed = subprocess.run([sys.executable, str(Path(__file__).with_name(script)), *paths, '--q-threshold', value], capture_output=True, text=True)
                self.assertEqual(completed.returncode, 2)
                self.assertIn('fixed 1% benchmark', completed.stderr)

    def test_fasta_mapping_labels(self):
        self.assertEqual(mapping_kind(self.row(ProteinName='KNNHHFK_target')), 'target')
        self.assertEqual(mapping_kind(self.row(ProteinName='LQREELR_p_target')), 'entrapment')
        self.assertEqual(mapping_kind(self.row(ProteinName='KNNHHFK_target;LQREELR_p_target')), 'ambiguous')


if __name__ == '__main__':
    unittest.main()
