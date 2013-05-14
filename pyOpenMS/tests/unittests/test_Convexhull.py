import pyopenms


def test_convex_hull_0():
    ch = pyopenms.ConvexHull2D()

def test_convex_hull_1():
    ch = pyopenms.ConvexHull2D()

    points = [ [1.0, 1.0], [2.0, 1.0], [1.0, 2.0], [2.0, 2.0]]
    for p in points:
        ch.addPoint(p)

    assert ch.encloses([1.5, 1.5])
    assert not ch.encloses([.5, 1.5])

    hp = ch.getHullPoints()
    # order may change, so we use sort here:
    assert  sorted(hp.tolist()) == sorted(points)

    ch.setHullPoints(hp)
    hp = ch.getHullPoints()
    # order may change, so we use sorted here:
    assert  sorted(hp.tolist()) == sorted(points)

    ch.addPoints(hp+1.0)
    ch.addPoints(hp)
    hp = ch.getHullPoints()
    assert hp.shape == (6, 2)

    bb = ch.getBoundingBox()
    assert bb.minPosition() == [1.0, 1.0]
    assert bb.maxPosition() == [3.0, 3.0]
