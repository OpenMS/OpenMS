import unittest
import os

import pyopenms

class TestDIAScoring(unittest.TestCase):

    def setUp(self):
        pass

    def test_spectrum(self):

          intensity = [
            100, 100, 100, 100,
            100, 100, 100]
          mz = [
                #// four of the naked b/y ions 
                #// as well as one of the modified b and y ions ion each
                350.17164, #// b
                421.20875, #// b
                421.20875 + 79.9657, #// b + P
                547.26291, #// y
                646.33133, #// y
                809.39466 + 79.9657 #// y + P
          ]

          spectrum = pyopenms.OSSpectrum()
          spectrum.setMZArray(mz)
          spectrum.setIntensityArray(intensity)

          diascoring = pyopenms.DIAScoring()
          # diascoring.set_dia_parameters(0.05, False, 30, 50, 4, 4) // here we use a large enough window so that none of our peaks falls out
          p_dia = diascoring.getDefaults();
          p_dia.setValue("dia_extraction_window", 0.05);
          p_dia.setValue("dia_extraction_unit", "Th");
          p_dia.setValue("dia_centroided", "false");
          p_dia.setValue("dia_byseries_intensity_min", 30.0);
          p_dia.setValue("dia_byseries_ppm_diff", 50.0);
          p_dia.setValue("dia_nr_isotopes", 4);
          p_dia.setValue("dia_nr_charges", 4);
          diascoring.setParameters(p_dia);

          a = pyopenms.AASequence.fromString(b"SYVAWDR", True)

          bseries_score = 0.0
          yseries_score = 0.0
          charge = 1
          bseries_score, yseries_score = diascoring.dia_by_ion_score(spectrum, a, charge, bseries_score, yseries_score)

          self.assertAlmostEqual(bseries_score, 2.0)
          self.assertAlmostEqual(yseries_score, 2.0)

          # // now add a modification to the sequence
          a.setModification(1, b"Phospho" ) #; // modify the Y
          bseries_score = 0
          yseries_score = 0
          bseries_score, yseries_score = diascoring.dia_by_ion_score(spectrum, a, 1, bseries_score, yseries_score) 

          self.assertAlmostEqual (bseries_score, 1.0)
          self.assertAlmostEqual (yseries_score, 3.0)

if __name__ == '__main__':
    unittest.main()
