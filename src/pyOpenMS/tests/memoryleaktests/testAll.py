from __future__ import print_function

import sys
import unittest
import os
import copy
import time
import contextlib
import pyopenms
from   pyopenms.sysinfo import free_mem
import numpy as np


def show_mem(label):

    p = free_mem()
    p /= 1024.0 * 1024
    print((label+" ").ljust(50, "."), ": %8.2f MB" % p)
    sys.stdout.flush()


@contextlib.contextmanager
def MemTester(name):
        mem_at_start = free_mem()
        print("\n")
        show_mem("start test '%s' with" % name)
        yield
        missing = mem_at_start - free_mem()
        show_mem("end with")
        print("\n")
        assert missing <= 0.1* mem_at_start, "possible mem leak: %s at start, missing %s" % (mem_at_start, missing)



if True or int(os.environ.get("WITH_MEMLEAK_TESTS", 0)):

    class TestAll(unittest.TestCase):

        def setUp(self):
            self.mem_at_start = free_mem()

            print("\n")
            show_mem("AT THE BEGINNING ")
            print("\n")

            dirname = os.path.dirname(os.path.abspath(__file__))
            self.testfile = os.path.join(dirname, "../test.mzXML").encode()

        def tearDown(self):

            time.sleep(3)
            print("\n")
            show_mem("AT THE END ")
            print("\n")
            missing = self.mem_at_start - free_mem()
            assert missing <= 0.1* self.mem_at_start, "possible mem leak: %s at start, missing %s" % (self.mem_at_start, missing)

        def testAll(self):

            with MemTester("specs from experiment"):
                self.run_extractSpetraFromMSExperiment()

            with MemTester("copy experiment"):
                self.run_MSExperiment_copy()

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

            basestr = 200000*b" "
            li = []
            for i in range(1000):
                if (i+1)%100 == 0:
                    show_mem("%4d runs" % i)
                dv = pyopenms.DataValue(basestr)
                dv = pyopenms.DataValue(basestr)
                li.append(dv)
            del li

        def run_string_conversions2(self):

            basestr = 200000*b" "
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

            basestr = 10000*b" "
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
                spec.set_peaks((data[:,0], data[:,1]))
                spec.set_peaks((data[:,0], data[:,1]))
                spec.set_peaks((data[:,0], data[:,1]))
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
                spec.set_peaks((data[:,0], data[:,1]))
                spec.set_peaks((data[:,0], data[:,1]))
                spec.set_peaks((data[:,0], data[:,1]))
                spec.set_peaks(spec.get_peaks())
                li.append(spec)

            for spec in li:
                del spec
            del data

        def run_extractSpetraFromMSExperiment(self):
            p = pyopenms.FileHandler()
            e = pyopenms.MSExperiment()
            p.loadExperiment(self.testfile, e)
            show_mem("data loaded")

            li = []
            print("please be patient :")
            for k in range(5):
                sys.stdout.flush()
                li.append([ e[i] for i in range(e.size()) ])
                li.append([ e[i] for i in range(e.size()) ])
                print((20*k+20), "%")

            print("\n")
            show_mem("spectra list generated")
            del li
            show_mem("spectra list deleted")
            del p
            del e

        def run_MSExperiment_copy(self):
            p = pyopenms.FileHandler()
            e1 = pyopenms.MSExperiment()
            p.loadExperiment(self.testfile, e1)
            show_mem("data loaded")

            specs = list(e1)
            for s in specs:
                for _ in range(10):
                    e1.addSpectrum(s)

            li = []
            print("please be patient :")
            N = 5
            for k in range(N):
                sys.stdout.flush()
                for __ in range(400):
                    e2 = copy.copy(e1)
                    li.append(e2)
                    e1 = copy.copy(e2)
                    li.append(e1)
                print(int(100.0*(k+1)/N), "%")

            print("\n")
            show_mem("experiment list generated")
            del li
            del e1
            del e2
            del p
            show_mem("experiment list deleted")

        def run_fileformats_io(self):
            p = pyopenms.FileHandler()
            e = pyopenms.MSExperiment()

            p.loadExperiment(self.testfile, e)
            show_mem("after load mzXML")

            ct = pyopenms.ChromatogramTools()
            ct.convertChromatogramsToSpectra(e)
            p.storeExperiment(self.testfile, e)
            show_mem("after store mzXML")

            p.loadExperiment(self.testfile, e)
            show_mem("after load mzXML")

            p = pyopenms.FileHandler()
            ct.convertSpectraToChromatograms(e, True, False)
            p.storeExperiment("../test.mzML".encode(), e)
            show_mem("after store mzML")
            p.loadExperiment("../test.mzML".encode(), e)
            show_mem("after load mzML")

            p = pyopenms.FileHandler()
            ct.convertChromatogramsToSpectra(e)
            p.storeExperiment("../test.mzData".encode(), e)
            show_mem("after store mzData")
            p.loadExperiment("../test.mzData".encode(), e)
            show_mem("after load mzData")

            del e
            del p
            del ct


if __name__ == "__main__":
    unittest.main()
