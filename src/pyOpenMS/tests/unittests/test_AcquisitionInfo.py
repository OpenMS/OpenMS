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

		self.assertEqual(exp[0].getAcquisitionInfo().size(), 1)
		self.assertNotEqual(exp[0].getAcquisitionInfo().size(), 0)
		self.assertEqual(exp[0].getAcquisitionInfo()[0].isMetaEmpty(), True)

if __name__ == '__main__':
	unittest.main()
