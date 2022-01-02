import unittest
import os

import pyopenms
class TestAcquisitionInfo(unittest.TestCase):

	def setUp(self):
		dirname = os.path.dirname(os.path.abspath(__file__))
		self.filename_mzml = os.path.join(dirname, "test.mzML").encode()

	def test_acquisitioninfomemberaccess(self):
		exp = pyopenms.MSExperiment()
		pyopenms.MzMLFile().load(self.filename_mzml, exp)

		# Basically test that the output is non-zero (e.g. the data is
		# correctly relayed to python)
		# The functionality is not tested here!

		# starting point
		self.assertEqual(exp[0].getAcquisitionInfo().size(), 1)
		self.assertNotEqual(exp[0].getAcquisitionInfo().size(), 0)

		# metainfo
		exp[0].getAcquisitionInfo().size()  # is 1
		self.assertEqual(exp[0].getAcquisitionInfo()[0].isMetaEmpty(), True)  # is True

		spectra = exp.getSpectra()
		aqis = spectra[0].getAcquisitionInfo()
		aqi = aqis[0] # get a copy
		aqi.setMetaValue('key', 420) # modify it
		aqis[0] = aqi # and set entry
		spectra[0].setAcquisitionInfo(aqis)
		exp.setSpectra(spectra)
		self.assertEqual(exp[0].getAcquisitionInfo()[0].getMetaValue('key'), 420)  # should be 420

		acin = pyopenms.Acquisition()
		acin.setMetaValue('key', 42)
		self.assertEqual(acin.getMetaValue('key'), 42)  # is 42
		self.assertEqual(acin.isMetaEmpty(), False)  # is False

		# list/vector assignment
		magicnumber = 3
		neac = pyopenms.AcquisitionInfo()
		for i in range(0,magicnumber):
			neac.push_back(acin)

		self.assertEqual(neac.size(), magicnumber)  # is magicnumber

		# iteration
		for i in neac:
			self.assertEqual(i.isMetaEmpty(), False)  # always is False
		
		# accession already tested in 2nd section
		tmp = exp.getSpectra()
		tmp[0].setAcquisitionInfo(neac)
		exp.setSpectra(tmp)

		self.assertEqual(exp[0].getAcquisitionInfo().size(), magicnumber)  # should be magicnumber

		for i in exp[0].getAcquisitionInfo():
			self.assertEqual(i.isMetaEmpty(), False)  # should always be False

		# resize
		neac.resize(0)
		self.assertEqual(neac.size(), 0)

if __name__ == '__main__':
	unittest.main()
