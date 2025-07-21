#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## $Maintainer: $
## $Authors: Hendrik Weisser, Hannes Roest, Stephan Aiche,
##           Jeremi Maciejewski $
## ----------------------------------------------------------------------------
import pyopenms
from functools import wraps

def report(f):
    @wraps(f)
    def wrapper(*a, **kw):
        print("run ", f.__name__)
        f(*a, **kw)
    return wrapper

@report
def testProteaseDigestion():
    """
    @tests: ProteaseDigestion
     ProteaseDigestion.__init__
     ProteaseDigestion.setEnzyme()
     ProteaseDigestion.digest()
     ProteaseDigestion.peptideCount()
     ProteaseDigestion.isValidProduct()
     ProteaseDigestion.getMissedCleavages()
     ProteaseDigestion.setMissedCleavages()
    """

    ff = pyopenms.ProteaseDigestion()

    assert pyopenms.ProteaseDigestion().setEnzyme is not None
    ff.setEnzyme("Trypsin/P")
    assert ff.getEnzymeName() == "Trypsin/P"

    assert pyopenms.ProteaseDigestion().digest is not None
    seq = pyopenms.AASequence.fromString("MHARLVP")
    output = []
    ff.digest(seq, output)
    # MHAR, LVP
    assert len(output) == 2
    ff.setSpecificity(1) # Semi-specific
    output = []
    ff.digest(seq, output)
    # MHAR, LVP, HAR, LV, AR, L, R, MHA, VP, MH, P, M
    assert len(output) == 12

    assert pyopenms.ProteaseDigestion().peptideCount is not None
    ff.setSpecificity(2)
    assert ff.peptideCount() == 2

    assert pyopenms.ProteaseDigestion().isValidProduct is not None
    assert ff.isValidProduct(seq, 0, 4) == True
    assert ff.isValidProduct(seq, 3, 4) == False
    ff.setSpecificity(1)
    assert ff.isValidProduct(seq, 1, 3) == True

    assert pyopenms.ProteaseDigestion().getMissedCleavages is not None
    assert pyopenms.ProteaseDigestion().setMissedCleavages is not None

    ff.setMissedCleavages(5)
    assert ff.getMissedCleavages() == 5