import unittest
import os

import pyopenms

class TestBilinearInterpolation(unittest.TestCase):

    def test_init(self):
        bilip = pyopenms.BilinearInterpolation()

    #TODO: def test_getData(self):
    #TODO: def test_setData(self):
    
    def test_setMapping_0(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setMapping_0(1.0, 2.0, 3.0, 8.0)
        
        self.assertAlmostEqual(bilip.getScale_0(), 3)
        self.assertAlmostEqual(bilip.getOffset_0(), -1)
        self.assertAlmostEqual(bilip.getScale_1(), 1)
        self.assertAlmostEqual(bilip.getOffset_1(), 0)
        
        bilip2 = pyopenms.BilinearInterpolation()
        bilip2.setMapping_0(3.0, 1.0, 2.0)
        
        self.assertAlmostEqual(bilip2.getScale_0(), 3)
        self.assertAlmostEqual(bilip2.getOffset_0(), -1)
        self.assertAlmostEqual(bilip2.getScale_1(), 1)
        self.assertAlmostEqual(bilip2.getOffset_1(), 0)      
        
    def test_setOffset_0(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setOffset_0(987)
        self.assertAlmostEqual(bilip.getOffset_0(), 987)
        
    def test_setScale_0(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setScale_0(987)
        self.assertAlmostEqual(bilip.getScale_0(), 987)
        
    def test_getInsideReferencePoint_0(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setMapping_0(1.0, 4.0, 3.0, 8.0)
        
        self.assertAlmostEqual(bilip.getInsideReferencePoint_0(), 1)
        self.assertAlmostEqual(bilip.getOutsideReferencePoint_0(), 4)
        self.assertAlmostEqual(bilip.getInsideReferencePoint_1(), 0)
        self.assertAlmostEqual(bilip.getOutsideReferencePoint_1(), 0)
        
    def test_index2key_0(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setMapping_0(3.0, 1.0, 2.0)
        self.assertAlmostEqual(bilip.index2key_0(0), -1)
        self.assertAlmostEqual(bilip.index2key_1(0), 0)
        
    def test_key2index_0(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setMapping_0(3.0, 1.0, 2.0)
        self.assertAlmostEqual(bilip.key2index_0(-1), 0)
        self.assertAlmostEqual(bilip.key2index_1(0), 0)
        
    def test_supportMax_0(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setMapping_0(3.0, 1.0, 2.0)
        bilip.setMapping_1(5.0, 3.0, 4.0)

        bilip.getData().resize(2, 3, float())
        self.assertAlmostEqual(bilip.index2key_0(0), -1)
        self.assertAlmostEqual(bilip.index2key_0(1), 2)
        self.assertAlmostEqual(bilip.supportMin_0(), -4)
        self.assertAlmostEqual(bilip.supportMax_0(), 5)

        self.assertAlmostEqual(bilip.index2key_1(0), -11)
        self.assertAlmostEqual(bilip.index2key_1(2), -1)
        self.assertAlmostEqual(bilip.supportMin_1(), -16)
        self.assertAlmostEqual(bilip.supportMax_1(), 4)
    
    def test_setMapping_1(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setMapping_1(1.0, 2.0, 3.0, 8.0)
        
        self.assertAlmostEqual(bilip.getScale_1(), 3)
        self.assertAlmostEqual(bilip.getOffset_1(), -1)
        self.assertAlmostEqual(bilip.getScale_0(), 1)
        self.assertAlmostEqual(bilip.getOffset_0(), 0)
    
        bilip2 = pyopenms.BilinearInterpolation()
        bilip2.setMapping_1(3.0, 1.0, 2.0)
        
        self.assertAlmostEqual(bilip2.getScale_1(), 3)
        self.assertAlmostEqual(bilip2.getOffset_1(), -1)
        self.assertAlmostEqual(bilip2.getScale_0(), 1)
        self.assertAlmostEqual(bilip2.getOffset_0(), 0)
        
    def test_setOffset_1(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setOffset_1(987)
        self.assertAlmostEqual(bilip.getOffset_1(), 987)

    def test_setScale_1(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setScale_1(987)
        self.assertAlmostEqual(bilip.getScale_1(), 987)
        
    def test_getInsideReferencePoint_1(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setMapping_1(1.0, 4.0, 3.0, 8.0)
        
        self.assertAlmostEqual(bilip.getInsideReferencePoint_1(), 1)
        self.assertAlmostEqual(bilip.getOutsideReferencePoint_1(), 4)
        self.assertAlmostEqual(bilip.getInsideReferencePoint_0(), 0)
        self.assertAlmostEqual(bilip.getOutsideReferencePoint_0(), 0)
        
    def test_index2key_1(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setMapping_1(3.0, 1.0, 2.0)
        self.assertAlmostEqual(bilip.index2key_1(0), -1)
        self.assertAlmostEqual(bilip.index2key_0(0), 0)

    def test_key2index_1(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setMapping_1(3.0, 1.0, 2.0)
        self.assertAlmostEqual(bilip.key2index_1(-1), 0)
        self.assertAlmostEqual(bilip.key2index_0(0), 0)

    def test_supportMax_1(self):
        bilip = pyopenms.BilinearInterpolation()
        bilip.setMapping_1(3.0, 1.0, 2.0)
        bilip.setMapping_0(5.0, 3.0, 4.0)

        bilip.getData().resize(3, 2, float())
        self.assertAlmostEqual(bilip.index2key_1(0), -1)
        self.assertAlmostEqual(bilip.index2key_1(1), 2)
        self.assertAlmostEqual(bilip.supportMin_1(), -4)
        self.assertAlmostEqual(bilip.supportMax_1(), 5)

        self.assertAlmostEqual(bilip.index2key_0(0), -11)
        self.assertAlmostEqual(bilip.index2key_0(2), -1)
        self.assertAlmostEqual(bilip.supportMin_0(), -16)
        self.assertAlmostEqual(bilip.supportMax_0(), 4)
        
    def test_empty(self):
        bilip = pyopenms.BilinearInterpolation()
        self.assertTrue(bilip.empty())

        bilip.getData().resize(1, 2, float())
        self.assertFalse(bilip.empty())

        bilip.getData().resize(0, 0, float())
        self.assertTrue(bilip.empty())

        bilip.getData().resize(1, 2, float())
        self.assertFalse(bilip.empty())

        bilip.getData().resize(1, 0, float())
        self.assertTrue(bilip.empty())

        bilip.getData().resize(1, 2, float())
        self.assertFalse(bilip.empty())

        bilip.getData().resize(0, 0, float())
        self.assertTrue(bilip.empty())

        bilip.getData().resize(2, 2, float())
        self.assertFalse(bilip.empty())

        bilip.getData().clear()
        self.assertTrue(bilip.empty())
        
    #TODO: test_addValue(self):
    #TODO: test_value(self):

if __name__ == '__main__':
    unittest.main()
    