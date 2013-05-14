import sys
import unittest
import time
import contextlib
import pyopenms
from   pyopenms.sysinfo import free_mem
import numpy as np


def show_mem(label):

    p = free_mem()
    p /= 1024.0 * 1024
    print (label+" ").ljust(50, "."), ": %8.2f MB" % p
    sys.stdout.flush()


@contextlib.contextmanager
def MemTester(name):
        mem_at_start = free_mem()
        print
        show_mem("start test '%s' with" % name)
        yield
        missing = mem_at_start - free_mem()
        show_mem("end with")
        print
        assert missing < 0.1* mem_at_start, "possible mem leak"


import os

if int(os.environ.get("WITH_MEMLEAK_TESTS", 0)):

    class TestAll(unittest.TestCase):

        def setUp(self):
            self.mem_at_start = free_mem()

            print
            show_mem("AT THE BEGINNING ")
            print

        def tearDown(self):

            time.sleep(3)
            print
            show_mem("AT THE END ")
            print
            missing = self.mem_at_start - free_mem()
            assert missing < 0.1* self.mem_at_start, "possible mem leak"

        def testAll(self):

            with MemTester("specs from experiment"):
                self.run_extractSpetraFromMSExperiment()

            with MemTester("string_conversions1"):
                self.run_string_conversions1()

            with MemTester("string_conversions2"):
                self.run_string_conversions2()

            with MemTester("string_lists"):
                self.run_string_lists()

            with MemTester("list_conversions"):
                self.run_list_conversions()

            with MemTester("set_spec_peaks"):
                self.set_spec_peaks()

            with MemTester("set_spec_peaks2"):
                self.set_spec_peaks2()


            with MemTester("test io"):
                self.run_fileformats_io()


        def run_string_conversions1(self):

            basestr = 200000*" "
            li = []
            for i in range(1000):
                if (i+1)%100 == 0:
                    show_mem("%4d runs" % i)
                dv = pyopenms.DataValue(basestr)
                dv = pyopenms.DataValue(basestr)
                li.append(dv)
            del li

        def run_string_conversions2(self):

            basestr = 200000*" "
            li = []
            for i in range(1000):
                if (i+1)%100 == 0:
                    show_mem("%4d runs" % i)
                sf = pyopenms.SourceFile()
                sf.setNameOfFile(basestr)
                sf.setNameOfFile(basestr)
                li.append(sf)
            del li

        def run_string_lists(self):

            basestr = 10000*" "
            li = []
            for i in range(100):
                if (i+1)%100 == 0:
                    show_mem("%4d runs" % i)
                sl = pyopenms.DataValue([basestr]*30)
                sl = pyopenms.DataValue([basestr]*30)
                li.append(sl)
                del sl
            del li

        def run_list_conversions(self):

            pc = pyopenms.Precursor()
            allpcs = 500*[pc]
            li = []
            for i in range(500):
                if (i+1)%100 == 0:
                    show_mem("%4d runs" % i)
                spec = pyopenms.MSSpectrum()
                spec.setPrecursors(allpcs)
                spec.setPrecursors(allpcs)
                li.append(spec)
                del spec
            del li

        def set_spec_peaks(self):

            data = np.zeros((10000,2), dtype=np.float32)
            li = []
            for i in range(1000):
                if (i+1)%100 == 0:
                    show_mem("%4d specs processed" % i)
                spec = pyopenms.MSSpectrum()
                spec.set_peaks(data)
                spec.set_peaks(data)
                spec.set_peaks(data)
                li.append(spec)

            for spec in li:
                del spec
            del data

        def set_spec_peaks2(self):

            data = np.zeros((10000,2), dtype=np.float32)
            li = []
            for i in range(1000):
                if (i+1)%100 == 0:
                    show_mem("%4d specs processed" % i)
                spec = pyopenms.MSSpectrum()
                spec.set_peaks(data)
                spec.set_peaks(data)
                spec.set_peaks(data)
                spec.set_peaks(spec.get_peaks())
                li.append(spec)

            for spec in li:
                del spec
            del data

        def run_extractSpetraFromMSExperiment(self):
                p = pyopenms.FileHandler()
                e = pyopenms.MSExperiment()
                p.loadExperiment("test.mzXML", e)
                show_mem("data loaded")

                li = []
                print "please be patient :",
                for k in range(5):
                    sys.stdout.flush()
                    li.append([ e[i] for i in range(e.size()) ])
                    li.append([ e[i] for i in range(e.size()) ])
                    print (20*k+20), "%",

                print
                show_mem("spectra list generated")
                del li
                show_mem("spectra list deleted")
                del p
                del e

        def run_fileformats_io(self):
            p = pyopenms.FileHandler()
            e = pyopenms.MSExperiment()

            p.loadExperiment("test.mzXML", e)
            show_mem("after load mzXML")

            ct = pyopenms.ChromatogramTools()
            ct.convertChromatogramsToSpectra(e)
            p.storeExperiment("test.mzXML", e)
            show_mem("after store mzXML")

            p.loadExperiment("test.mzXML", e)
            show_mem("after load mzXML")

            p = pyopenms.FileHandler()
            ct.convertSpectraToChromatograms(e, True)
            p.storeExperiment("test.mzML", e)
            show_mem("after store mzML")
            p.loadExperiment("test.mzML", e)
            show_mem("after load mzML")

            p = pyopenms.FileHandler()
            ct.convertChromatogramsToSpectra(e)
            p.storeExperiment("test.mzData", e)
            show_mem("after store mzData")
            p.loadExperiment("test.mzData", e)
            show_mem("after load mzData")

            del e
            del p
            del ct


if __name__ == "__main__":
    unittest.main()
