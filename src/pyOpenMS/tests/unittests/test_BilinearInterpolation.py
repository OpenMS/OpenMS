import unittest
import os

import pyopenms

class TestBilinearInterpolation(unittest.TestCase):

    def test_init(self):
        bilip = pyopenms.BilinearInterpolation()


if __name__ == '__main__':
    unittest.main()
    