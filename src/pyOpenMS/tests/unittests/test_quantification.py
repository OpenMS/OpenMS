#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## Quantification tests extracted from test000.py
## Part of Issue #8567: Split test000.py into modular test files
## ----------------------------------------------------------------------------

from __future__ import print_function
import pyopenms
import os

# Import shared test helper functions
from test_helpers import report, _testParam


@report
def testItraqConstants():
  """
  @tests: testItraqConstants
  """
  constants = pyopenms.ItraqConstants()

  assert pyopenms.ITRAQ_TYPES.FOURPLEX is not None
  assert pyopenms.ITRAQ_TYPES.EIGHTPLEX is not None
  assert pyopenms.ITRAQ_TYPES.TMT_SIXPLEX is not None

  assert constants.getIsotopeMatrixAsStringList is not None
  assert constants.updateIsotopeMatrixFromStringList is not None
  assert constants.translateIsotopeMatrix is not None


@report
def testPeptideAndProteinQuant():
  """
  @tests: PeptideAndProteinQuant
   PeptideAndProteinQuant.__init__
  """
  ff = pyopenms.PeptideAndProteinQuant()
  p = ff.getDefaults()
  _testParam(p)

  assert pyopenms.PeptideAndProteinQuant().quantifyPeptides is not None
  assert pyopenms.PeptideAndProteinQuant().quantifyProteins is not None


@report
def testGNPSExport():
  """
  @tests: GNPSExport
   IonIdentityMolecularNetworking.annotateConsensusMap
   IonIdentityMolecularNetworking.writeSupplementaryPairTable
   GNPSQuantificationFile.store
   GNPSMetaValueFile.store
  """
  cm = pyopenms.ConsensusMap()
    
  for mz, rt, ion, linked_groups in [(222.08, 62.0, "[M+H]+", ["1","2","3"]),
                    (244.08, 62.0, "[M+Na]+", ["1","2"]),
                    (204.08, 62.0, "[M-H-O]+", ["3"]),
                    (294.1, 62.0, "[M+H]+", ["4","5"])]:
    f = pyopenms.ConsensusFeature()
    f.setMZ(mz)
    f.setRT(rt)
    f.setCharge(1)
    f.setQuality(2.0)
    f.setMetaValue("best ion", ion)
    f.setMetaValue("LinkedGroups", linked_groups)
    cm.push_back(f)
  cm.setUniqueIds()

  pyopenms.IonIdentityMolecularNetworking.annotateConsensusMap(cm)

  pyopenms.IonIdentityMolecularNetworking.writeSupplementaryPairTable(cm, "SupplementaryPairsTable.csv")

  with open("SupplementaryPairsTable.csv", "r") as f:
    assert f.read() == """ID1,ID2,EdgeType,Score,Annotation
1,2,MS1 annotation,1,[M+H]+ [M+Na]+ dm/z=22.0
1,3,MS1 annotation,1,[M+H]+ [M-H-O]+ dm/z=18.0
"""
  os.remove("SupplementaryPairsTable.csv")

  pyopenms.GNPSQuantificationFile().store(cm, "FeatureQuantificationTable.txt")
  with open("FeatureQuantificationTable.txt", "r") as f:
    assert f.read() == """#MAP	id	filename	label	size
#CONSENSUS	rt_cf	mz_cf	intensity_cf	charge_cf	width_cf	quality_cf	row ID	best ion	partners	annotation network number
CONSENSUS	62.0	222.080000000000013	0.0	1	0.0	2.0	1	[M+H]+	2;3	1
CONSENSUS	62.0	244.080000000000013	0.0	1	0.0	2.0	2	[M+Na]+	1	1
CONSENSUS	62.0	204.080000000000013	0.0	1	0.0	2.0	3	[M-H-O]+	1	1
CONSENSUS	62.0	294.100000000000023	0.0	1	0.0	2.0	4	[M+H]+		2
"""
  os.remove("FeatureQuantificationTable.txt")

  # add mandatory file descriptions
  file_descriptions = {}
  for i, filename in enumerate(["1.mzML", "2.mzML"]):
    file_description = pyopenms.ColumnHeader()
    file_description.filename = filename
    file_descriptions[i] = file_description
  cm.setColumnHeaders(file_descriptions)

  pyopenms.GNPSMetaValueFile().store(cm, "MetaValueTable.tsv")
  with open("MetaValueTable.tsv", "r") as f:
    assert f.read() == """	filename	ATTRIBUTE_MAPID
0	1.mzML	MAP0
1	2.mzML	MAP1
"""
  os.remove("MetaValueTable.tsv")


if __name__ == "__main__":
  # Run all tests in this module
  import sys
  module = sys.modules[__name__]
    
  test_functions = [getattr(module, name) for name in dir(module) 
           if name.startswith('test') and callable(getattr(module, name))]
    
  print(f"Running {len(test_functions)} quantification tests...")
  for test_func in test_functions:
    try:
      test_func()
      print(f"✅ {test_func.__name__}")
    except Exception as e:
      print(f"❌ {test_func.__name__}: {e}")
