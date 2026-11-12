import unittest
import numpy as np
import pandas as pd
from helpers import winners,percentile,auc,membership_category
class Helpers(unittest.TestCase):
    def test_ties(self):
        np.testing.assert_allclose(percentile([0,1,2],[0,1,1,2]),[.125,.5,.875])
        self.assertEqual(auc([1,1],[1,1]),.5)
    def test_winner(self):
        p=pd.DataFrame({'precursor_id':[1,1,2],'sample_id':[9,2,3]})
        self.assertEqual(winners(p,[1,1,0]).tolist(),[1,2])
    def test_membership(self):
        self.assertEqual(membership_category([[0,0,0],[1,0,0],[0,1,1],[1,0,1]]).tolist(),['never_admitted','admitted_then_lost','retained_r2','retained_r2'])
if __name__=='__main__': unittest.main()
