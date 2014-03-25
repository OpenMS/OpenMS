import pyopenms
import os.path

from collections_ import Counter


def test0():
    fh = pyopenms.MzXMLFile()
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "test2.mzXML")

    class Consumer(object):

        def __init__(self):
            self.speclevels = []
            self.rts = []

        def consumeSpectrum(self, spec):
            self.speclevels.append(spec.getMSLevel())
            self.rts.append(spec.getRT())

        def consumeChromatogram(self, chromo):
            raise Exception("should never be called as we have no chromoatograms in example file")

        def setExpectedSize(self, num_specs, num_chromo):
            assert num_specs == 5, num_specs
            assert num_chromo == 0, num_chromo

        def setExperimentalSettings(self, exp):
            assert isinstance(exp, pyopenms.ExperimentalSettings)

    consumer = Consumer()
    fh.transform(path, consumer)
    cc = Counter(consumer.speclevels)
    assert set(cc.keys()) == set([1, 2])
    assert cc[1] == 2
    assert cc[2] == 3
    assert abs(min(consumer.rts) - 4200.76) < 0.01
    assert abs(max(consumer.rts) - 4202.03) < 0.01
