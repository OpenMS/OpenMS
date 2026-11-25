"""Synthetic tests only, never access baseline or scheduler."""
from pathlib import Path
import tempfile
import unittest

from compare_resources import artifacts, compare, evaluate, resources


class ComparisonTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name) / 'work'
        self.base = Path(self.temp.name) / 'base'
        self.accounts = {}
        for i in range(15):
            raw, seed = i // 5 + 1, i % 5
            name = f'raw{raw:02d}_s{seed}'
            for root in [self.work, self.base]:
                directory = root / 'runs' / name
                directory.mkdir(parents=True)
                (directory / 'exit_code.txt').write_text('0\n')
                for file in artifacts(raw):
                    (directory / file).write_bytes(b'identical\n')
            self.accounts[f'123_{i}'] = {'JobID': f'123_{i}', 'State': 'COMPLETED', 'ExitCode': '0:0'}

    def result(self):
        return evaluate(self.work, self.base, self.accounts, '123')

    def test_all15_match(self):
        result = self.result()
        self.assertEqual(result['denominator'], 15)
        self.assertEqual(result['passed'], 15)
        self.assertEqual([(r['raw'], r['seed']) for r in result['runs']],
                         [(raw, seed) for raw in range(1, 4) for seed in range(5)])

    def test_absent_all15_retained(self):
        result = evaluate(self.work / 'absent', self.base, {}, '123')
        self.assertEqual((result['denominator'], len(result['runs']), result['passed']), (15, 15, 0))

    def test_missing_each_output_fails(self):
        directory = self.work / 'runs/raw01_s0'
        for file in artifacts(1):
            with self.subTest(file=file):
                path = directory / file
                path.unlink()
                self.assertEqual(self.result()['passed'], 14)
                path.write_bytes(b'identical\n')

    def test_different_each_output_fails(self):
        directory = self.work / 'runs/raw01_s0'
        for file in artifacts(1):
            with self.subTest(file=file):
                path = directory / file
                path.write_bytes(b'different\n')
                self.assertEqual(self.result()['passed'], 14)
                path.write_bytes(b'identical\n')

    def test_failed_and_missing_actual_exit(self):
        path = self.work / 'runs/raw01_s0/exit_code.txt'
        path.write_text('9\n')
        self.assertEqual(self.result()['passed'], 14)
        path.unlink()
        (path.parent / 'complete').touch()
        self.assertEqual(self.result()['passed'], 14)

    def test_accounting_required(self):
        for state, code in [('FAILED', '1:0'), ('TIMEOUT', '0:0'), ('COMPLETED', '0:9')]:
            with self.subTest(state=state, code=code):
                self.accounts['123_0'].update(State=state, ExitCode=code)
                self.assertEqual(self.result()['passed'], 14)
        del self.accounts['123_0']
        self.assertEqual(self.result()['passed'], 14)

    def test_baseline_exit_required(self):
        (self.base / 'runs/raw01_s0/exit_code.txt').unlink()
        self.assertEqual(self.result()['passed'], 14)

    def test_bootstrap_uses_same_strict_comparison(self):
        directory = self.work / 'runs/raw01_s0'
        result = compare(directory, self.base / 'runs/raw01_s0', 1, '123_0', self.accounts)
        self.assertEqual(result['status'], 'PASS')

    def test_all15_failed_exits_retained(self):
        for path in (self.work / 'runs').glob('*/exit_code.txt'):
            path.write_text('1\n')
        result = self.result()
        self.assertEqual((len(result['runs']), result['denominator'], result['passed']), (15, 15, 0))

    def test_all15_different_retained(self):
        for path in (self.work / 'runs').glob('*/precursors.tsv'):
            path.write_bytes(b'different')
        result = self.result()
        self.assertEqual((len(result['runs']), result['denominator'], result['passed']), (15, 15, 0))

    def test_resources(self):
        directory = self.work / 'runs/raw01_s0'
        (directory / 'time.txt').write_text('Elapsed (wall clock) time (h:mm:ss or m:ss): 1:02:03.5\nMaximum resident set size (kbytes): 1234\n')
        self.assertEqual(resources(directory), {'wall_seconds': 3723.5, 'max_rss_kib': 1234})
        self.assertEqual(resources(self.work), {'wall_seconds': None, 'max_rss_kib': None})


class ShellContractTest(unittest.TestCase):
    directory = Path(__file__).parent

    def script(self, name):
        return (self.directory / name).read_text()

    def test_full_preflight_only_once_at_submit(self):
        self.assertEqual(self.script('submit.sh').count('sha256sum -c preflight.sha256'), 1)
        self.assertNotIn('preflight.sha256', self.script('run.sbatch'))

    def test_primary_uses_small_runtime_and_cache_checks(self):
        script = self.script('run.sbatch')
        self.assertIn('sha256sum -c runtime.sha256', script)
        self.assertIn('sha256sum -c library.sha256', script)
        self.assertNotIn('sacct ', script)
        self.assertNotIn('.pqp', script)
        self.assertIn('-lib "$WORK/library/search6.opendia.parquet"', script)

    def test_library_manifest_is_cache_only_and_checked_before_scheduling(self):
        script = self.script('submit.sh')
        self.assertIn('library/search6.opendia.parquet', script)
        self.assertIn('library/scoring12.opendia.parquet', script)
        self.assertNotIn('.pqp', script)
        self.assertLess(script.index('sha256sum -c library.sha256'), script.index('sbatch --parsable'))

    def test_runtime_manifest_required_before_scheduling(self):
        script = self.script('submit.sh')
        self.assertLess(script.index('test -s runtime.sha256'), script.index('sbatch --parsable'))
        self.assertLess(script.index('sha256sum -c runtime.sha256'), script.index('sbatch --parsable'))

    def test_final_accounting_remains_mandatory(self):
        script = self.script('evaluate.sbatch')
        self.assertIn('sacct -j "$(cat primary_job_id.txt)"', script)
        self.assertIn('"$B" > results/comparison.json', script)
        self.assertIn('"$C" > results/comparison_compact.json', script)
        self.assertIn('test "$primary_rc" = 0', script)
        self.assertIn('--dependency="afterany:$primary"', self.script('submit.sh'))


if __name__ == '__main__':
    unittest.main()
