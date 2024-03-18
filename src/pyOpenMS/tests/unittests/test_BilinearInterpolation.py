import unittest
import math
import random

import pyopenms


class TestBilinearInterpolation(unittest.TestCase):

    def test_BilinearInterpolation(self):
        mat = pyopenms.MatrixDouble(2, 3, 0.0)
        mat.setValue(0, 0, 17)
        mat.setValue(0, 1, 18.9)
        mat.setValue(0, 2, 20.333)
        mat.setValue(1, 0, -0.1)
        mat.setValue(1, 1, -0.13)
        mat.setValue(1, 2, -0.001)

        bilip = pyopenms.BilinearInterpolation()
        bilip.setData(mat)
        bilip.setMapping_0(13.0, 230.0, 14.0, 250.0)
        bilip.setMapping_1(15.0, 2100.0, 17.0, 2900.0)

        bilip_2 = pyopenms.BilinearInterpolation()
        bilip_2 = bilip

        self.assertAlmostEqual(bilip_2.getScale_0(), 20)
        self.assertAlmostEqual(bilip_2.getScale_1(), 400)
        self.assertAlmostEqual(bilip_2.getOffset_0(), -30)
        self.assertAlmostEqual(bilip_2.getOffset_1(), -3900)
        self.assertAlmostEqual(bilip_2.getInsideReferencePoint_0(), 13)
        self.assertAlmostEqual(bilip_2.getOutsideReferencePoint_0(), 230)
        self.assertAlmostEqual(bilip_2.getInsideReferencePoint_1(), 15)
        self.assertAlmostEqual(bilip_2.getOutsideReferencePoint_1(), 2100)

        for i in range(mat.rows()):
            for j in range(mat.cols()):
                self.assertAlmostEqual(bilip.getData().getValue(i, j),
                                       mat.getValue(i, j))

    # same test for getData and setData
    def test_getData_setData(self):
        bilip = pyopenms.BilinearInterpolation()

        tmp = pyopenms.MatrixDouble(2, 3, 0.0)

        tmp.setValue(1, 2, 10012)
        tmp.setValue(0, 0, 10000)
        tmp.setValue(1, 0, 10010)
        bilip.setData(tmp)

        bilip_2 = pyopenms.BilinearInterpolation(bilip)
        self.assertAlmostEqual(bilip_2.getData().getValue(1, 2), 10012)
        self.assertAlmostEqual(bilip_2.getData().getValue(0, 0), 10000)
        self.assertAlmostEqual(bilip_2.getData().getValue(1, 0), 10010)

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

        tmp = pyopenms.MatrixDouble(2, 3, 0.0)

        bilip.setData(tmp)
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

        tmp = pyopenms.MatrixDouble(3, 2, 0.0)
        bilip.setData(tmp)
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

        tmp = pyopenms.MatrixDouble(1, 2, 0.0)
        bilip.setData(tmp)
        self.assertFalse(bilip.empty())

        tmp = pyopenms.MatrixDouble(0, 0, 0.0)
        bilip.setData(tmp)
        self.assertTrue(bilip.empty())

        tmp = pyopenms.MatrixDouble(1, 2, 0.0)
        bilip.setData(tmp)
        self.assertFalse(bilip.empty())

        tmp = pyopenms.MatrixDouble(1, 0, 0.0)
        bilip.setData(tmp)
        self.assertTrue(bilip.empty())

        tmp = pyopenms.MatrixDouble(1, 2, 0.0)
        bilip.setData(tmp)
        self.assertFalse(bilip.empty())

        tmp = pyopenms.MatrixDouble(0, 0, 0.0)
        bilip.setData(tmp)
        self.assertTrue(bilip.empty())

        tmp = pyopenms.MatrixDouble(2, 2, 0.0)
        bilip.setData(tmp)
        self.assertFalse(bilip.empty())

        tmp = bilip.getData()
        tmp = pyopenms.MatrixDouble(0, 0, 0.0)
        bilip.setData(tmp)
        self.assertTrue(bilip.empty())

    def test_addValue(self):
        for i in range(-50, 101):
            p = i / 10.0

            for j in range(-50, 101):
                q = j / 10.0

                bilip_small = pyopenms.BilinearInterpolation()
                tmp = pyopenms.MatrixDouble(5, 5, 0.0)
                bilip_small.setData(tmp)
                bilip_small.setMapping_0(0.0, 0.0, 5.0, 5.0)
                bilip_small.setMapping_1(0.0, 0.0, 5.0, 5.0)
                bilip_small.addValue(p, q, 100)

                bilip_big = pyopenms.BilinearInterpolation()
                tmp = pyopenms.MatrixDouble(15, 15, 0.0)
                bilip_big.setData(tmp)
                bilip_big.setMapping_0(5.0, 0.0, 10.0, 5.0)
                bilip_big.setMapping_1(5.0, 0.0, 10.0, 5.0)
                bilip_big.addValue(p, q, 100)

                big_submatrix = pyopenms.MatrixDouble(5, 5, 0.0)

                for m in range(5):
                    for n in range(5):
                        big_submatrix.setValue(m, n,
                                               bilip_big.getData().getValue(m+5,
                                                                            n+5))

                for m in range(bilip_small.getData().rows()):
                    for n in range(bilip_small.getData().cols()):
                        self.assertAlmostEqual(bilip_small.getData().getValue(m, n),
                                               big_submatrix.getValue(m, n))

    def test_value(self):
        bilip_small = pyopenms.BilinearInterpolation()
        bilip_big = pyopenms.BilinearInterpolation()

        tmp_small = pyopenms.MatrixDouble(5, 5, 0.0)
        tmp_big = pyopenms.MatrixDouble(15, 15, 0.0)

        for i in range(5):
            for j in range(5):
                num = int(math.floor(random.random() * 100))
                tmp_small.setValue(i, j, num)
                tmp_big.setValue(i + 5, j + 5, num)

        bilip_small.setData(tmp_small)
        bilip_big.setData(tmp_big)

        bilip_small.setMapping_0(0.0, 0.0, 5.0, 5.0)
        bilip_small.setMapping_1(0.0, 0.0, 5.0, 5.0)
        bilip_big.setMapping_0(5.0, 0.0, 10.0, 5.0)
        bilip_big.setMapping_1(5.0, 0.0, 10.0, 5.0)

        for i in range(-50, 101):
            p = i / 10.0
            for j in range(-50, 101):
                q = j / 10.0
                self.assertAlmostEqual(bilip_small.value(p, q),
                                       bilip_big.value(p, q))

if __name__ == '__main__':
    unittest.main()
