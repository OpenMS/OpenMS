import pyopenms
import numpy

def _test_container(c, i):
    # did raise segfault because of missing index check in __getitem__
    # should be resolved by implementing a proper iterator:
    s = c()
    assert list(s) == []
    try:
        s[0]
    except IndexError:
        pass
    else:
        assert False, "no exception risen"

    s.push_back(i())
    (item,) = s # test iterator

    s[0]  # test getattr

    try:
        s[1]
    except IndexError:
        pass
    else:
        assert False, "no exception risen"

def testContainers():
    _test_container(pyopenms.MSSpectrum, pyopenms.Peak1D)
    _test_container(pyopenms.MSExperiment, pyopenms.MSSpectrum)
    _test_container(pyopenms.FeatureMap, pyopenms.Feature)


def testParamEntry():
    # as ParamEntry::isValid takes "String &" as input argument, which
    # can not be implemened by a Python string, here no automatic
    # conversion from a basestring should happen:
    p = pyopenms.ParamEntry()
    message = pyopenms.String()
    assert p.isValid(message)
    assert message.c_str() == ""


def testConvexHull2D():
    h = pyopenms.ConvexHull2D()

    points = numpy.arange(10.0, dtype=numpy.float32).reshape(-1, 2)

    h.setHullPoints(points)
    p = h.getHullPoints()
    assert numpy.linalg.norm(p - points) < 1e-6

    h.expandToBoundingBox()
    hullp = h.getHullPoints()
    assert set(hullp[:,0]) == set((0.0, 8.0))
    assert set(hullp[:,1]) == set((1.0, 9.0))

    box = h.getBoundingBox()
    assert box.minPosition() == [ 0.0, 1.0]
    assert box.maxPosition() == [ 8.0, 9.0]

    h.addPoint([-10.0, -10.0])
    box = h.getBoundingBox()
    assert box.minPosition() == [ -10.0, -10.0]
    assert box.maxPosition() == [ 8.0, 9.0]
    h.expandToBoundingBox()
    hullp = h.getHullPoints()
    assert set(hullp[:,0]) == set((-10.0, 8.0))
    assert set(hullp[:,1]) == set((-10.0, 9.0))
    box = h.getBoundingBox()
    assert box.minPosition() == [ -10.0, -10.0]
    assert box.maxPosition() == [ 8.0, 9.0]

