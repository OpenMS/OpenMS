"""Integration tests for all imzML load modes (mirrors ImzMLFile_all_modes_test.cpp)."""

import os
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from test_data_paths import get_class_test_data_dir

import pyopenms


def _repo_imzml_paths():
    try:
        data_dir = get_class_test_data_dir()
    except RuntimeError as e:
        raise unittest.SkipTest(str(e)) from e
    continuous = os.path.join(data_dir, "ImzMLFile_1_Example_Continuous.imzML")
    processed = os.path.join(data_dir, "ImzMLFile_2_Example_Processed.imzML")
    return continuous, processed


def _find_spectrum_index_by_imzml_coord(exp, x, y):
    for i in range(exp.getNrSpectra()):
        spec = exp.getSpectrum(i)
        if not spec.metaValueExists("imzml:x") or not spec.metaValueExists("imzml:y"):
            continue
        if int(spec.getMetaValue("imzml:x")) == x and int(spec.getMetaValue("imzml:y")) == y:
            return i
    raise AssertionError(f"No spectrum at imzML pixel ({x}, {y})")


class TestImzMLAllModes(unittest.TestCase):
    CONTINUOUS_SPECTRA = 9
    GRID = 3

    @classmethod
    def setUpClass(cls):
        cls.continuous_path, cls.processed_path = _repo_imzml_paths()
        if not os.path.isfile(cls.continuous_path):
            raise unittest.SkipTest(f"imzML test data not found: {cls.continuous_path}")

    def test_mode1_full_load_msexperiment(self):
        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.continuous_path, exp)

        self.assertEqual(exp.getNrSpectra(), self.CONTINUOUS_SPECTRA)
        self.assertGreater(exp.getSpectrum(0).size(), 0)
        self.assertTrue(exp.metaValueExists("imzml:imaging_mode"))
        self.assertEqual(str(exp.getMetaValue("imzml:imaging_mode")), "continuous")
        self.assertEqual(int(exp.getMetaValue("imzml:max_count_x")), self.GRID)
        self.assertEqual(int(exp.getMetaValue("imzml:max_count_y")), self.GRID)
        self.assertEqual(int(exp.getSpectrum(0).getMetaValue("imzml:x")), 1)
        self.assertEqual(int(exp.getSpectrum(0).getMetaValue("imzml:y")), 1)

    def test_mode2_file_handler_and_file_types(self):
        self.assertEqual(
            pyopenms.FileHandler.getType(self.continuous_path),
            pyopenms.FileType.IMZML,
        )
        self.assertEqual(
            pyopenms.FileHandler.getTypeByContent(self.continuous_path),
            pyopenms.FileType.IMZML,
        )
        exp = pyopenms.MSExperiment()
        pyopenms.FileHandler().loadExperiment(self.continuous_path, exp)
        self.assertEqual(exp.getNrSpectra(), self.CONTINUOUS_SPECTRA)
        self.assertGreater(exp.getSpectrum(0).size(), 0)

    def test_mode3_streaming_consumer(self):
        class Consumer:
            def __init__(self):
                self.count = 0
                self.first_size = 0

            def setExpectedSize(self, num_specs, num_chromo):
                pass

            def setExperimentalSettings(self, exp):
                assert isinstance(exp, pyopenms.ExperimentalSettings)

            def consumeChromatogram(self, chromo):
                pass

            def consumeSpectrum(self, spec):
                if self.count == 0:
                    self.first_size = spec.size()
                self.count += 1

        consumer = Consumer()
        pyopenms.ImzMLFile().load(self.continuous_path, consumer)
        self.assertEqual(consumer.count, self.CONTINUOUS_SPECTRA)
        self.assertGreater(consumer.first_size, 0)

    def test_mode4_load_spectra_index(self):
        meta, index = pyopenms.ImzMLFile().loadSpectraIndex(self.continuous_path)

        self.assertEqual(len(index), self.CONTINUOUS_SPECTRA)
        self.assertEqual(meta.imaging_mode, "continuous")
        self.assertEqual(meta.max_count_x, self.GRID)
        self.assertEqual(meta.max_count_y, self.GRID)
        if index:
            self.assertEqual(index[0].x, 1)
            self.assertEqual(index[0].y, 1)
            self.assertGreater(index[0].mz_length, 0)

    def test_mode5_on_disc_random_access(self):
        od = pyopenms.OnDiscImzMLExperiment()
        try:
            od.open(self.continuous_path)

            self.assertTrue(od.isOpen())
            self.assertEqual(od.getNrSpectra(), self.CONTINUOUS_SPECTRA)
            self.assertEqual(od.gridWidth(), self.GRID)
            self.assertEqual(od.gridHeight(), self.GRID)
            self.assertEqual(od.getImzMLMeta().imaging_mode, "continuous")

            by_index = od.getSpectrum(0)
            self.assertGreater(by_index.size(), 0)
            e0 = od.getIndex(0)
            by_coord = od.getSpectrumAtCoord(e0.x, e0.y, e0.z)
            self.assertEqual(by_coord.size(), by_index.size())
        finally:
            od.close()

    def test_mode6_imaging_experiment_ram_lookup(self):
        imaging = pyopenms.MSImagingExperiment()
        pyopenms.ImzMLFile().load(self.continuous_path, imaging)

        self.assertEqual(imaging.getNumberOfSpectra(), self.CONTINUOUS_SPECTRA)
        self.assertEqual(imaging.getNumberOfPixels(), self.CONTINUOUS_SPECTRA)
        self.assertEqual(imaging.getGeometry().getWidth(), self.GRID)
        self.assertEqual(imaging.getGeometry().getHeight(), self.GRID)
        self.assertTrue(imaging.hasPixel(0, 0))
        self.assertTrue(imaging.hasPixel(2, 2))
        self.assertGreater(imaging.getSpectrum(0, 0).size(), 0)
        self.assertGreater(imaging.getSpectrum(2, 2).size(), 0)

    def test_mode7_build_imaging_geometry(self):
        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.continuous_path, exp)

        geom = pyopenms.MSImagingGeometry()
        pyopenms.ImzMLFile.buildImagingGeometry(exp, geom)
        self.assertEqual(geom.getNumberOfPixels(), self.CONTINUOUS_SPECTRA)
        self.assertEqual(geom.getSpectrumIndex(0, 0), 0)
        self.assertEqual(
            geom.getSpectrumIndex(2, 2),
            _find_spectrum_index_by_imzml_coord(exp, 3, 3),
        )

    def test_cross_mode_peak_count_consistency(self):
        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.continuous_path, exp)

        imaging = pyopenms.MSImagingExperiment()
        pyopenms.ImzMLFile().load(self.continuous_path, imaging)

        od = pyopenms.OnDiscImzMLExperiment()
        try:
            od.open(self.continuous_path)

            px, py = 2, 3
            idx = _find_spectrum_index_by_imzml_coord(exp, px, py)
            spec = exp.getSpectrum(idx)
            pz = int(spec.getMetaValue("imzml:z")) if spec.metaValueExists("imzml:z") else 1
            full_size = spec.size()
            ram_size = imaging.getSpectrum(px - 1, py - 1).size()
            disc_size = od.getSpectrumAtCoord(px, py, pz).size()

            self.assertGreater(full_size, 0)
            self.assertEqual(ram_size, full_size)
            self.assertEqual(disc_size, full_size)
        finally:
            od.close()

    def test_processed_encoding(self):
        if not os.path.isfile(self.processed_path):
            self.skipTest(f"processed imzML test data not found: {self.processed_path}")

        exp = pyopenms.MSExperiment()
        pyopenms.ImzMLFile().load(self.processed_path, exp)
        self.assertGreater(exp.getNrSpectra(), 0)
        self.assertGreater(exp.getSpectrum(0).size(), 0)
        self.assertEqual(str(exp.getMetaValue("imzml:imaging_mode")), "processed")

        od = pyopenms.OnDiscImzMLExperiment()
        try:
            od.open(self.processed_path)
            self.assertEqual(od.getImzMLMeta().imaging_mode, "processed")
            self.assertGreater(od.getSpectrum(0).size(), 0)
        finally:
            od.close()


if __name__ == "__main__":
    unittest.main()
