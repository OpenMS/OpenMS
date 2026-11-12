import unittest
import numpy as np
import pandas as pd
from mismatch_helpers import compare_rows
class SelectorHelpers(unittest.TestCase):
    def test_alignment_and_strict_distance(self):
        p=pd.DataFrame({'precursor_id':[1,2,1,2,3],'rt_apex':[0,2,10,13,4]})
        a,d,rt=compare_rows(p,[4,0,1],[3,4,2])
        np.testing.assert_array_equal(a,[4,2,3])
        np.testing.assert_array_equal(d,[False,True,True])
        np.testing.assert_array_equal(rt>10,[False,False,True])
    def test_duplicate_alternative_rejected(self):
        p=pd.DataFrame({'precursor_id':[1,1],'rt_apex':[0,2]})
        with self.assertRaises(AssertionError): compare_rows(p,[0],[0,1])
    def test_missing_alternative_rejected(self):
        p=pd.DataFrame({'precursor_id':[1,2],'rt_apex':[0,2]})
        with self.assertRaises(AssertionError): compare_rows(p,[0],[1])
if __name__=='__main__': unittest.main()
