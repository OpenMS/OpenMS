import unittest
import os

import pyopenms

class TestMSNumpressCoder(unittest.TestCase):

    def setUp(self):
        self.testData = [ 
            100.0,
            200.0,
            300.00005,
            400.00010,
        ]

    def test_encodeNP_SLOF(self):
        """
          String out;

          MSNumpressCoder::NumpressConfig config;
          config.np_compression = MSNumpressCoder::SLOF;
          config.estimate_fixed_point = true; // critical

          bool zlib_compression = false;
          MSNumpressCoder().encodeNP(in, out, zlib_compression, config);

          TEST_EQUAL(out.size(), 24)
          TEST_EQUAL(out, "QMVagAAAAAAZxX3ivPP8/w==")

        """
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.SLOF
        config.estimate_fixed_point = True

        out_ = pyopenms.String()
        coder.encodeNP(self.testData, out_, False, config)
        out = out_.c_str()

        self.assertEqual( len(out),  24)
        self.assertEqual( out,  b"QMVagAAAAAAZxX3ivPP8/w==")

    def test_decodeNP_SLOF(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.SLOF
        config.estimate_fixed_point = True

        inData = b"QMVagAAAAAAZxX3ivPP8/w=="
        out = []
        coder.decodeNP(inData, out, False, config)

        self.assertEqual( len(out),  4)
        for a,b in zip(self.testData, out):
            self.assertAlmostEqual( a, b, places=2)

    def test_encodeNP_PIC(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.PIC
        config.estimate_fixed_point = True

        out_ = pyopenms.String()
        coder.encodeNP(self.testData, out_, False, config)
        out = out_.c_str()

        self.assertEqual( len(out),  12)
        self.assertEqual( out, b"ZGaMXCFQkQ==")

    def test_decodeNP_PIC(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.PIC
        config.estimate_fixed_point = True

        inData = b"ZGaMXCFQkQ=="
        out = []
        coder.decodeNP(inData, out, False, config)

        self.assertEqual( len(out),  4)
        for a,b in zip(self.testData, out):
            self.assertAlmostEqual( a, b, places=2)

    def test_encodeNP_LINEAR(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.LINEAR
        config.estimate_fixed_point = True

        out_ = pyopenms.String()
        coder.encodeNP(self.testData, out_, False, config)
        out = out_.c_str()

        self.assertEqual( len(out),  28)
        self.assertEqual( out,  b"QWR64UAAAADo//8/0P//f1kSgA==")

    def test_decodeNP_LINEAR(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.LINEAR
        config.estimate_fixed_point = True

        inData = b"QWR64UAAAADo//8/0P//f1kSgA=="
        out = []
        coder.decodeNP(inData, out, False, config)

        self.assertEqual( len(out),  4)
        for a,b in zip(self.testData, out):
            self.assertAlmostEqual( a, b, places=7)

class TestMSNumpressCoderRaw(unittest.TestCase):

    def setUp(self):
        self.testData = [ 
            100.0,
            200.0,
            300.00005,
            400.00010,
        ]

    def test_encodeNP_SLOF(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.SLOF
        config.estimate_fixed_point = True

        out_ = pyopenms.String()
        coder.encodeNPRaw(self.testData, out_, config)
        out = out_.c_str()

        self.assertEqual( len(out),  16)
        self.assertEqual( out,  b'@\xc5Z\x80\x00\x00\x00\x00\x19\xc5}\xe2\xbc\xf3\xfc\xff' )

    def test_decodeNP_SLOF(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.SLOF
        config.estimate_fixed_point = True

        inData = b'@\xc5Z\x80\x00\x00\x00\x00\x19\xc5}\xe2\xbc\xf3\xfc\xff'
        out = []
        coder.decodeNPRaw(inData, out, config)

        self.assertEqual( len(out),  4)
        for a,b in zip(self.testData, out):
            self.assertAlmostEqual( a, b, places=2)

    def test_encodeNP_PIC(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.PIC
        config.estimate_fixed_point = True

        out_ = pyopenms.String()
        coder.encodeNPRaw(self.testData, out_, config)
        out = out_.c_str()

        self.assertEqual( len(out),  7)
        self.assertEqual( out, b'df\x8c\\!P\x91')

    def test_decodeNP_PIC(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.PIC
        config.estimate_fixed_point = True

        inData = b'df\x8c\\!P\x91'
        out = []
        coder.decodeNPRaw(inData, out, config)

        self.assertEqual( len(out),  4)
        for a,b in zip(self.testData, out):
            self.assertAlmostEqual( a, b, places=2)

    def test_encodeNP_LINEAR(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.LINEAR
        config.estimate_fixed_point = True

        out_ = pyopenms.String()
        coder.encodeNPRaw(self.testData, out_, config)
        out = out_.c_str()

        self.assertEqual( len(out),  19)
        self.assertEqual( out,  b'Adz\xe1@\x00\x00\x00\xe8\xff\xff?\xd0\xff\xff\x7fY\x12\x80')

    def test_decodeNP_LINEAR(self):
        coder = pyopenms.MSNumpressCoder()
        config = pyopenms.NumpressConfig()
        config.np_compression = pyopenms.MSNumpressCoder.NumpressCompression.LINEAR
        config.estimate_fixed_point = True

        inData = b'Adz\xe1@\x00\x00\x00\xe8\xff\xff?\xd0\xff\xff\x7fY\x12\x80'
        out = []
        coder.decodeNPRaw(inData, out, config)

        self.assertEqual( len(out),  4)
        for a,b in zip(self.testData, out):
            self.assertAlmostEqual( a, b, places=7)

if __name__ == '__main__':
    unittest.main()
