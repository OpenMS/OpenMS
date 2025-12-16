"""
Tests for DRange1 Python wrapper
"""

import unittest
import pyopenms


class TestDRange1(unittest.TestCase):
    """Test DRange1 class"""

    def test_default_constructor(self):
        """Test default constructor creates empty range"""
        r = pyopenms.DRange1()
        self.assertTrue(r.isEmpty())

    def test_copy_constructor(self):
        """Test copy constructor"""
        r1 = pyopenms.DRange1(10.0, 20.0)
        r2 = pyopenms.DRange1(r1)
        self.assertEqual(r1.minX(), r2.minX())
        self.assertEqual(r1.maxX(), r2.maxX())

    def test_two_float_constructor(self):
        """Test constructor with two floats"""
        r = pyopenms.DRange1(10.0, 20.0)
        self.assertFalse(r.isEmpty())
        self.assertEqual(r.minX(), 10.0)
        self.assertEqual(r.maxX(), 20.0)

    def test_two_int_constructor(self):
        """Test constructor accepts ints and converts to floats"""
        r = pyopenms.DRange1(5, 15)
        self.assertEqual(r.minX(), 5.0)
        self.assertEqual(r.maxX(), 15.0)

    def test_setMinX_setMaxX(self):
        """Test setMinX and setMaxX methods"""
        r = pyopenms.DRange1()
        r.setMinX(100.0)
        r.setMaxX(200.0)
        self.assertEqual(r.minX(), 100.0)
        self.assertEqual(r.maxX(), 200.0)

    def test_getMin_getMax(self):
        """Test getMin/getMax aliases"""
        r = pyopenms.DRange1(10.0, 20.0)
        self.assertEqual(r.getMin(), 10.0)
        self.assertEqual(r.getMax(), 20.0)

    def test_encloses(self):
        """Test encloses method with float values"""
        r = pyopenms.DRange1(10.0, 20.0)

        # Below min - not enclosed
        self.assertFalse(r.encloses(5.0))

        # At min - enclosed (inclusive)
        self.assertTrue(r.encloses(10.0))

        # Inside range - enclosed
        self.assertTrue(r.encloses(15.0))

        # At max - not enclosed (exclusive)
        self.assertFalse(r.encloses(20.0))

        # Above max - not enclosed
        self.assertFalse(r.encloses(25.0))

    def test_united(self):
        """Test united method"""
        r1 = pyopenms.DRange1(5.0, 15.0)
        r2 = pyopenms.DRange1(12.0, 25.0)
        r3 = r1.united(r2)
        self.assertEqual(r3.minX(), 5.0)
        self.assertEqual(r3.maxX(), 25.0)

    def test_isIntersected(self):
        """Test isIntersected method"""
        r1 = pyopenms.DRange1(5.0, 15.0)
        r2 = pyopenms.DRange1(12.0, 25.0)
        r3 = pyopenms.DRange1(20.0, 30.0)

        # Overlapping ranges
        self.assertTrue(r1.isIntersected(r2))

        # Non-overlapping ranges
        self.assertFalse(r1.isIntersected(r3))

    def test_repr_str(self):
        """Test __repr__ and __str__ methods"""
        r = pyopenms.DRange1(10.0, 20.0)
        self.assertIn("DRange1", repr(r))
        self.assertIn("10.0", repr(r))
        self.assertIn("20.0", repr(r))
        self.assertEqual(str(r), repr(r))


class TestPeakFileOptionsWithDRange1(unittest.TestCase):
    """Test PeakFileOptions methods that use DRange1"""

    def test_rt_range(self):
        """Test setRTRange/getRTRange/hasRTRange"""
        pfo = pyopenms.PeakFileOptions()

        self.assertFalse(pfo.hasRTRange())

        rt_range = pyopenms.DRange1(100.0, 500.0)
        pfo.setRTRange(rt_range)

        self.assertTrue(pfo.hasRTRange())

        rt_back = pfo.getRTRange()
        self.assertEqual(rt_back.minX(), 100.0)
        self.assertEqual(rt_back.maxX(), 500.0)

    def test_mz_range(self):
        """Test setMZRange/getMZRange/hasMZRange"""
        pfo = pyopenms.PeakFileOptions()

        self.assertFalse(pfo.hasMZRange())

        mz_range = pyopenms.DRange1(200.0, 1000.0)
        pfo.setMZRange(mz_range)

        self.assertTrue(pfo.hasMZRange())

        mz_back = pfo.getMZRange()
        self.assertEqual(mz_back.minX(), 200.0)
        self.assertEqual(mz_back.maxX(), 1000.0)

    def test_intensity_range(self):
        """Test setIntensityRange/getIntensityRange/hasIntensityRange"""
        pfo = pyopenms.PeakFileOptions()

        self.assertFalse(pfo.hasIntensityRange())

        int_range = pyopenms.DRange1(1000.0, 1e6)
        pfo.setIntensityRange(int_range)

        self.assertTrue(pfo.hasIntensityRange())

        int_back = pfo.getIntensityRange()
        self.assertEqual(int_back.minX(), 1000.0)
        self.assertEqual(int_back.maxX(), 1e6)


if __name__ == "__main__":
    unittest.main()
