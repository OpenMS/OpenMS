import os
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from test_data_paths import get_class_test_data_dir

import pyopenms


class TestImzMLFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            data_dir = get_class_test_data_dir()
        except RuntimeError as e:
            raise unittest.SkipTest(str(e)) from e
        cls.imzml_path = os.path.join(
            data_dir,
            "ImzMLFile_1_Example_Continuous.imzML",
        )
        cls.processed_path = os.path.join(
            data_dir,
            "ImzMLFile_2_Example_Processed.imzML",
        )
        if not os.path.isfile(cls.imzml_path):
            raise unittest.SkipTest(f"imzML test data not found: {cls.imzml_path}")

    @staticmethod
    def _make_pixel_spectrum(x, y, mz=100.0, intensity=1000.0):
        spec = pyopenms.MSSpectrum()
        peak = pyopenms.Peak1D()
        peak.setMZ(mz)
        peak.setIntensity(intensity)
        spec.push_back(peak)
        spec.setMetaValue("imzml:x", x)
        spec.setMetaValue("imzml:y", y)
        spec.setMetaValue("imzml:z", 1)
        return spec

    def test_imzml_file_load(self):
        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.imzml_path, exp)
        self.assertGreater(exp.getNrSpectra(), 0)
        self.assertGreater(exp.getSpectrum(0).size(), 0)

    def test_file_handler_load(self):
        exp = pyopenms.MSExperiment()
        pyopenms.FileHandler().loadExperiment(self.imzml_path, exp)
        self.assertGreater(exp.getNrSpectra(), 0)

    def test_file_types(self):
        t = pyopenms.FileHandler.getType(self.imzml_path)
        self.assertEqual(t, pyopenms.FileType.IMZML)

    def test_on_disc(self):
        od = pyopenms.OnDiscImzMLExperiment()
        try:
            od.open(self.imzml_path)
            self.assertTrue(od.isOpen())
            self.assertGreater(od.getNrSpectra(), 0)
            s = od.getSpectrum(0)
            self.assertGreater(s.size(), 0)
        finally:
            od.close()
            self.assertFalse(od.isOpen())

    def test_imaging_experiment_ram_lookup(self):
        imaging = pyopenms.MSImagingExperiment()
        pyopenms.ImzMLFile().load(self.imzml_path, imaging)
        self.assertEqual(imaging.getGeometry().getWidth(), 3)
        self.assertEqual(imaging.getGeometry().getHeight(), 3)
        self.assertTrue(imaging.hasPixel(0, 0))
        self.assertGreater(imaging.getSpectrum(0, 0).size(), 0)
        self.assertTrue(imaging.hasPixel(2, 2))
        self.assertGreater(imaging.getSpectrum(2, 2).size(), 0)

    def test_load_spectra_index(self):
        meta, index = pyopenms.ImzMLFile().loadSpectraIndex(self.imzml_path)
        self.assertEqual(meta.imaging_mode, "continuous")
        self.assertEqual(meta.max_count_x, 3)
        self.assertEqual(meta.max_count_y, 3)
        self.assertEqual(len(index), 9)
        self.assertEqual(index[0].x, 1)
        self.assertEqual(index[0].y, 1)
        self.assertGreater(index[0].mz_length, 0)

    def test_on_disc_index_and_meta(self):
        od = pyopenms.OnDiscImzMLExperiment()
        try:
            od.open(self.imzml_path)
            meta = od.getImzMLMeta()
            self.assertEqual(meta.imaging_mode, "continuous")
            entry = od.getIndex(0)
            self.assertEqual(entry.x, 1)
            self.assertEqual(entry.y, 1)
            self.assertGreater(entry.mz_length, 0)
        finally:
            od.close()

    def test_store_round_trip(self):
        import tempfile

        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.imzml_path, exp)
        with tempfile.NamedTemporaryFile(suffix=".imzML", delete=False) as tmp:
            out_path = tmp.name
        try:
            pyopenms.ImzMLFile().store(out_path, exp)
            reloaded = pyopenms.MSExperiment()
            pyopenms.ImzMLFile().load(out_path, reloaded)
            self.assertEqual(reloaded.getNrSpectra(), exp.getNrSpectra())
            self.assertGreater(reloaded.getSpectrum(0).size(), 0)
            self.assertEqual(str(reloaded.getMetaValue("imzml:imaging_mode")), "continuous")
            self.assertTrue(reloaded.metaValueExists("imzml:ibd_md5"))
            self.assertTrue(str(reloaded.getMetaValue("imzml:ibd_md5")))
        finally:
            os.remove(out_path)
            ibd_path = out_path[:-6] + ".ibd" if out_path.lower().endswith(".imzml") else out_path + ".ibd"
            if os.path.isfile(ibd_path):
                os.remove(ibd_path)

    def test_store_round_trip_processed(self):
        if not os.path.isfile(self.processed_path):
            self.skipTest(f"processed imzML test data not found: {self.processed_path}")

        import tempfile

        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.processed_path, exp)
        with tempfile.NamedTemporaryFile(suffix=".imzML", delete=False) as tmp:
            out_path = tmp.name
        try:
            pyopenms.ImzMLFile().store(out_path, exp)
            reloaded = pyopenms.MSExperiment()
            pyopenms.ImzMLFile().load(out_path, reloaded)
            self.assertEqual(reloaded.getNrSpectra(), exp.getNrSpectra())
            self.assertGreater(reloaded.getSpectrum(0).size(), 0)
            self.assertEqual(str(reloaded.getMetaValue("imzml:imaging_mode")), "processed")
        finally:
            os.remove(out_path)
            ibd_path = out_path[:-6] + ".ibd" if out_path.lower().endswith(".imzml") else out_path + ".ibd"
            if os.path.isfile(ibd_path):
                os.remove(ibd_path)

    def test_is_valid(self):
        import tempfile

        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.imzml_path, exp)
        with tempfile.NamedTemporaryFile(suffix=".imzML", delete=False) as tmp:
            out_path = tmp.name
        try:
            pyopenms.ImzMLFile().store(out_path, exp)
            ok, _msg = pyopenms.ImzMLFile().isValid(out_path)
            self.assertTrue(ok)
        finally:
            os.remove(out_path)
            ibd_path = out_path[:-6] + ".ibd" if out_path.lower().endswith(".imzml") else out_path + ".ibd"
            if os.path.isfile(ibd_path):
                os.remove(ibd_path)

    def test_store_rejects_missing_pixel_coordinates(self):
        import tempfile

        exp = pyopenms.MSExperiment()
        spec = pyopenms.MSSpectrum()
        peak = pyopenms.Peak1D()
        peak.setMZ(100.0)
        peak.setIntensity(1000.0)
        spec.push_back(peak)
        exp.addSpectrum(spec)

        with tempfile.NamedTemporaryFile(suffix=".imzML", delete=False) as tmp:
            out_path = tmp.name
        try:
            with self.assertRaises(Exception):
                pyopenms.ImzMLFile().store(out_path, exp)
        finally:
            if os.path.isfile(out_path):
                os.remove(out_path)

    def test_store_rejects_duplicate_pixel_coordinates(self):
        import tempfile

        exp = pyopenms.MSExperiment()
        exp.addSpectrum(self._make_pixel_spectrum(1, 1, 100.0, 1000.0))
        exp.addSpectrum(self._make_pixel_spectrum(1, 1, 101.0, 900.0))

        with tempfile.NamedTemporaryFile(suffix=".imzML", delete=False) as tmp:
            out_path = tmp.name
        try:
            with self.assertRaises(Exception):
                pyopenms.ImzMLFile().store(out_path, exp)
        finally:
            if os.path.isfile(out_path):
                os.remove(out_path)

    def test_store_rejects_incompatible_continuous_mode(self):
        import tempfile

        exp = pyopenms.MSExperiment()
        exp.setMetaValue("imzml:imaging_mode", "continuous")
        exp.addSpectrum(self._make_pixel_spectrum(1, 1, 100.0, 1000.0))
        exp.addSpectrum(self._make_pixel_spectrum(2, 1, 200.0, 800.0))

        with tempfile.NamedTemporaryFile(suffix=".imzML", delete=False) as tmp:
            out_path = tmp.name
        try:
            with self.assertRaises(Exception):
                pyopenms.ImzMLFile().store(out_path, exp)
        finally:
            if os.path.isfile(out_path):
                os.remove(out_path)

    def test_build_imaging_geometry_rejects_duplicate_pixels(self):
        exp = pyopenms.MSExperiment()
        exp.addSpectrum(self._make_pixel_spectrum(1, 1, 100.0, 1000.0))
        exp.addSpectrum(self._make_pixel_spectrum(1, 1, 101.0, 900.0))
        geom = pyopenms.MSImagingGeometry()
        with self.assertRaises(Exception):
            pyopenms.ImzMLFile.buildImagingGeometry(exp, geom)

    def test_store_applies_mz_range_filter(self):
        import tempfile

        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.imzml_path, exp)
        full_peaks = exp.getSpectrum(0).size()
        lo = exp.getSpectrum(0)[0].getMZ()
        hi = exp.getSpectrum(0)[full_peaks - 1].getMZ()
        mid = (lo + hi) / 2.0

        f = pyopenms.ImzMLFile()
        opts = pyopenms.PeakFileOptions()
        opts.setMZRange(pyopenms.DRange1(lo, mid))
        f.setOptions(opts)

        with tempfile.NamedTemporaryFile(suffix=".imzML", delete=False) as tmp:
            out_path = tmp.name
        try:
            f.store(out_path, exp)
            reloaded = pyopenms.MSExperiment()
            pyopenms.ImzMLFile().load(out_path, reloaded)
            self.assertLess(reloaded.getSpectrum(0).size(), full_peaks)
            self.assertLessEqual(reloaded.getSpectrum(0)[reloaded.getSpectrum(0).size() - 1].getMZ(), mid)
        finally:
            os.remove(out_path)
            ibd_path = out_path[:-6] + ".ibd" if out_path.lower().endswith(".imzml") else out_path + ".ibd"
            if os.path.isfile(ibd_path):
                os.remove(ibd_path)

    def test_store_binary_precision_options(self):
        import tempfile

        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.imzml_path, exp)

        f = pyopenms.ImzMLFile()
        opts = pyopenms.PeakFileOptions()
        opts.setMz32Bit(False)
        opts.setIntensity32Bit(False)
        f.setOptions(opts)

        with tempfile.NamedTemporaryFile(suffix=".imzML", delete=False) as tmp:
            out_path = tmp.name
        try:
            f.store(out_path, exp)
            _meta, index = pyopenms.ImzMLFile().loadSpectraIndex(out_path)
            self.assertEqual(index[0].mz_type, pyopenms.ImzMLDataType.FLOAT64)
            self.assertEqual(index[0].int_type, pyopenms.ImzMLDataType.FLOAT64)
        finally:
            os.remove(out_path)
            ibd_path = out_path[:-6] + ".ibd" if out_path.lower().endswith(".imzml") else out_path + ".ibd"
            if os.path.isfile(ibd_path):
                os.remove(ibd_path)

    def test_store_metadata_round_trip(self):
        import tempfile

        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.imzml_path, exp)
        exp.setMetaValue("imzml:scan_pattern", "top down")
        exp.setMetaValue("imzml:polarity", "positive")
        inst = exp.getInstrument()
        inst.setModel("Test MSI Instrument")
        exp.setInstrument(inst)
        if exp.getNrSpectra() > 0:
            spec = exp.getSpectrum(0)
            spec.setRT(12.34)
            spec.setMSLevel(1)

        with tempfile.NamedTemporaryFile(suffix=".imzML", delete=False) as tmp:
            out_path = tmp.name
        try:
            pyopenms.ImzMLFile().store(out_path, exp)
            reloaded = pyopenms.MSExperiment()
            pyopenms.ImzMLFile().load(out_path, reloaded)
            self.assertEqual(str(reloaded.getMetaValue("imzml:scan_pattern")), "top down")
            self.assertEqual(str(reloaded.getMetaValue("imzml:polarity")), "positive")
            self.assertTrue(reloaded.metaValueExists("imzml:ibd_sha1"))
            self.assertTrue(str(reloaded.getMetaValue("imzml:ibd_sha1")))
            if reloaded.getNrSpectra() > 0:
                self.assertAlmostEqual(reloaded.getSpectrum(0).getRT(), 12.34, places=4)
                self.assertEqual(reloaded.getSpectrum(0).getMSLevel(), 1)
        finally:
            os.remove(out_path)
            ibd_path = out_path[:-6] + ".ibd" if out_path.lower().endswith(".imzml") else out_path + ".ibd"
            if os.path.isfile(ibd_path):
                os.remove(ibd_path)


if __name__ == "__main__":
    unittest.main()
