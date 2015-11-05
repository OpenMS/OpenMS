// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer:	David Wojnar $
// $Authors: David Wojnar, Petra Gutenbrunner $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>


///////////////////////////
#include <OpenMS/ANALYSIS/ID/AScore.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class AScoreTest : public AScore
{
  public:
    void computeSiteDeterminingIonsTest_(std::vector<RichPeakSpectrum> & th_spectra, ProbablePhosphoSites & candidates, std::vector<RichPeakSpectrum> & site_determining_ions) const
    {
      return computeSiteDeterminingIons_(th_spectra, candidates, site_determining_ions);
    }

    std::vector<Size> getSitesTest_(AASequence & without_phospho) const {
      return getSites_(without_phospho);
    }

    std::vector<std::vector<Size> > computePermutationsTest_(std::vector<Size> sites, Int n_phosphorylation_events) const
    {
      return computePermutations_(sites, n_phosphorylation_events);
    }

    Size numberOfMatchedIonsTest_(const RichPeakSpectrum & th, const RichPeakSpectrum & windows, Size depth, double fragment_mass_tolerance, bool fragment_mass_tolerance_ppm = false) const
    {
      return numberOfMatchedIons_(th, windows, depth, fragment_mass_tolerance, fragment_mass_tolerance_ppm);
    }

    //double peptideScoreTest_(const std::vector<double> & scores) const;

    void determineHighestScoringPermutationsTest_(const std::vector<std::vector<double> > & peptide_site_scores, std::vector<ProbablePhosphoSites> & sites, const std::vector<std::vector<Size> > & permutations, std::multimap<double, Size>& ranking) const
    {
      return determineHighestScoringPermutations_(peptide_site_scores, sites, permutations, ranking);
    }
    
    double computeCumulativeScoreTest_(Size N, Size n, double p) const
    {
      return computeCumulativeScore_(N, n, p);
    }  
    
    //Size numberOfPhosphoEventsTest_(const String sequence) const;
    
    AASequence removePhosphositesFromSequenceTest_(const String sequence) const
    {
      return removePhosphositesFromSequence_(sequence);
    }
    
    std::vector<RichPeakSpectrum> createTheoreticalSpectraTest_(const std::vector<std::vector<Size> > & permutations, const AASequence & seq_without_phospho) const
    {
      return createTheoreticalSpectra_(permutations, seq_without_phospho);
    }
    
    std::vector<RichPeakSpectrum> createTheoreticalSpectraTest_(const AASequence & seq_without_phospho) const
    {
      return createTheoreticalSpectra_(seq_without_phospho);
    }
    
    std::vector<RichPeakSpectrum> peakPickingPerWindowsInSpectrumTest_(RichPeakSpectrum & real_spectrum) const
    {
      return peakPickingPerWindowsInSpectrum_(real_spectrum);
    }
    
    //std::vector<std::vector<double> > calculatePermutationPeptideScoresTest_(std::vector<RichPeakSpectrum> & th_spectra, const std::vector<RichPeakSpectrum> & windows_top10, double fragment_mass_tolerance, bool fragment_mass_unit_ppm) const;
    
    std::multimap<double, Size> rankWeightedPermutationPeptideScoresTest_(const std::vector<std::vector<double> > & peptide_site_scores) const
    {
      return rankWeightedPermutationPeptideScores_(peptide_site_scores);
    }
};

///////////////////////////
///////////////////////////

START_TEST(AScore, "$Id$")

//=============================================================================
// peak data see Beausoleil et al. Figure 3
//=============================================================================
  
//b3
RichPeak1D peak6;
peak6.getPosition()[0] = 303.1;
peak6.setIntensity(1.0f);

//y3
RichPeak1D peak;
peak.getPosition()[0] = 411.1;
peak.setIntensity(1.0f);

//y4
RichPeak1D peak2;
peak2.getPosition()[0] = 539.2;
peak2.setIntensity(1.0f);

//y5
RichPeak1D peak3;
peak3.getPosition()[0] = 668.3;
peak3.setIntensity(1.0f);

//b9
RichPeak1D peak4;
peak4.getPosition()[0] = 960.5;
peak4.setIntensity(1.0f);

//b10
RichPeak1D peak5;
peak5.getPosition()[0] = 1088.5;
peak5.setIntensity(1.0f);

//=============================================================================
// create spectrum
//=============================================================================
MSSpectrum<RichPeak1D > tmp;  
tmp.push_back(peak6);
tmp.push_back(peak);
tmp.push_back(peak2);
tmp.push_back(peak3);
tmp.push_back(peak4);
tmp.push_back(peak5);

//=============================================================================
  AASequence seq_without_phospho = AASequence::fromString("QSSVTQVTEQSPK");
//=============================================================================

//=============================================================================
// create permutations based on sequence QSSVTQVTEQSPK
//=============================================================================
std::vector<std::vector<Size> > permutations;
std::vector<Size> perm;

perm.clear();
perm.push_back(1);
permutations.push_back(perm);

perm.clear();
perm.push_back(2);
permutations.push_back(perm);

perm.clear();
perm.push_back(4);
permutations.push_back(perm);

perm.clear();
perm.push_back(7);
permutations.push_back(perm);

perm.clear();
perm.push_back(10);
permutations.push_back(perm);
//=============================================================================

AScore* ptr = 0;
AScore* nullPointer = 0;
START_SECTION(AScore())
{
  ptr = new AScore();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~AScore())
{
  delete ptr;
}
END_SECTION

AScoreTest* ptr_test = new AScoreTest();

START_SECTION(double computeCumulativeScoreTest_(Size N, Size n, double p))
{
  Size n = 5;
  Size N = 1;
  double p = 0.1;
  TEST_PRECONDITION_VIOLATED(ptr_test->computeCumulativeScoreTest_(N,n,p));

  n = 1;
  double score = ptr_test->computeCumulativeScoreTest_(N,n,p);
  TEST_REAL_SIMILAR(score,0.1);
  N = 3;
  score = ptr_test->computeCumulativeScoreTest_(N,n,p);
  TEST_REAL_SIMILAR(score,0.271);
}
END_SECTION

START_SECTION(determineHighestScoringPermutationsTest_(const std::vector<std::vector<double> > & peptide_site_scores, std::vector<ProbablePhosphoSites> & sites, const std::vector<std::vector<Size> > & permutations))
{
  std::multimap<double, Size> ranking;
  std::vector< std::vector<double> > peptide_site_scores_1;
  std::vector< std::vector<double> > peptide_site_scores_2;
  std::vector< std::vector<double> > peptide_site_scores_3;
  peptide_site_scores_1.resize(4);
  peptide_site_scores_2.resize(4);
  peptide_site_scores_3.resize(4);
  vector<double> temp;
  temp.resize(10);
  for(Size i = 0; i < 10; ++i)
  {
    temp[i] = 0.1;
  }
  peptide_site_scores_1[0] = temp;
  peptide_site_scores_2[3] = temp;
  peptide_site_scores_3[0] = temp;
  temp.clear();
  temp.resize(10);
  for(Size i = 0; i < 10; ++i)
  {
    temp[i] = 0.2;
  }
  peptide_site_scores_1[1] = temp;
  peptide_site_scores_2[0] = temp;
  peptide_site_scores_3[3] = temp;
  temp.clear();
  temp.resize(10);
  for(Size i = 0; i < 10; ++i)
  {
    temp[i] = 0.3;
  }
  peptide_site_scores_1[2] = temp;
  peptide_site_scores_2[1] = temp;
  peptide_site_scores_3[2] = temp;
  temp.clear();
  temp.resize(10);
  for(Size i = 0; i < 10; ++i)
  {
    temp[i] = 0.4;
  }
  peptide_site_scores_1[3] = temp;
  peptide_site_scores_2[2] = temp;
  peptide_site_scores_3[1] = temp;


  vector<vector<Size> > permutations;
  vector<Size> per;
  per.push_back(1);
  per.push_back(3);
  per.push_back(5);
  permutations.push_back(per);
  per.clear();
  per.push_back(3);
  per.push_back(5);
  per.push_back(6);
  permutations.push_back(per);
  per.clear();
  per.push_back(1);
  per.push_back(3);
  per.push_back(6);
  permutations.push_back(per);
  per.clear();
  per.push_back(1);
  per.push_back(5);
  per.push_back(6);
  permutations.push_back(per);


  vector<ProbablePhosphoSites> sites;
  ranking = ptr_test->rankWeightedPermutationPeptideScoresTest_(peptide_site_scores_1);
  ptr_test->determineHighestScoringPermutationsTest_(peptide_site_scores_1,sites,permutations,ranking);
  TEST_EQUAL(sites.size(),3)
  TEST_EQUAL(sites[0].seq_1, 3);
  TEST_EQUAL(sites[0].seq_2,1);
  TEST_EQUAL(sites[0].second, 3);
  TEST_EQUAL(sites[0].first,1);
  TEST_EQUAL(sites[0].peak_depth, 1)
    TEST_EQUAL(sites[1].first,5);
  TEST_EQUAL(sites[1].second,3);
  TEST_EQUAL(sites[1].seq_1, 3);
  TEST_EQUAL(sites[1].seq_2,2);
  TEST_EQUAL(sites[1].peak_depth, 1)
    TEST_EQUAL(sites[2].first,6);
  TEST_EQUAL(sites[2].second,3);
    TEST_EQUAL(sites[2].seq_1, 3);
  TEST_EQUAL(sites[2].seq_2,0);
  TEST_EQUAL(sites[2].peak_depth, 1)

  ranking = ptr_test->rankWeightedPermutationPeptideScoresTest_(peptide_site_scores_3);
  ptr_test->determineHighestScoringPermutationsTest_(peptide_site_scores_3,sites,permutations,ranking);
  TEST_EQUAL(sites.size(),3)
  TEST_EQUAL(sites[0].seq_1, 1);
  TEST_EQUAL(sites[0].seq_2,3);
  TEST_EQUAL(sites[0].second,1 );
  TEST_EQUAL(sites[0].first,3);
  TEST_EQUAL(sites[0].peak_depth, 1)
    TEST_EQUAL(sites[1].first,5);
  TEST_EQUAL(sites[1].second,1);
  TEST_EQUAL(sites[1].seq_1, 1);
  TEST_EQUAL(sites[1].seq_2,2);
  TEST_EQUAL(sites[1].peak_depth, 1)
    TEST_EQUAL(sites[2].first,6);
  TEST_EQUAL(sites[2].second,1);
    TEST_EQUAL(sites[2].seq_1, 1);
  TEST_EQUAL(sites[2].seq_2,0);
  TEST_EQUAL(sites[2].peak_depth, 1)

  ranking = ptr_test->rankWeightedPermutationPeptideScoresTest_(peptide_site_scores_2);
  ptr_test->determineHighestScoringPermutationsTest_(peptide_site_scores_2,sites,permutations,ranking);
  TEST_EQUAL(sites.size(),3)
  TEST_EQUAL(sites[0].seq_1, 2);
  TEST_EQUAL(sites[0].seq_2,1);
  TEST_EQUAL(sites[0].second,5 );
  TEST_EQUAL(sites[0].first,1);
  TEST_EQUAL(sites[0].peak_depth, 1)
    TEST_EQUAL(sites[1].first,3);
  TEST_EQUAL(sites[1].second,5);
  TEST_EQUAL(sites[1].seq_1, 2);
  TEST_EQUAL(sites[1].seq_2,3);
  TEST_EQUAL(sites[1].peak_depth, 1)
    TEST_EQUAL(sites[2].first,6);
  TEST_EQUAL(sites[2].second,5);
    TEST_EQUAL(sites[2].seq_1, 2);
  TEST_EQUAL(sites[2].seq_2,0);
  TEST_EQUAL(sites[2].peak_depth, 1)

  peptide_site_scores_1.clear();
  temp.clear();
  temp.resize(10);
  temp[0] = 55;
  temp[1] = 60;
  temp[2]= 75;
  temp[3] = 100;
  temp[4] = 90;
  temp[5] = 120;
  temp[6]  =125;
  temp[7] = 120;
  temp[8] = 100;
  temp[9] = 90;
  peptide_site_scores_1.push_back(temp);
  temp.clear();
  temp.resize(10);
  temp[0] = 40;
  temp[1] = 50;
  temp[2] = 53;
  temp[3] = 60;
  temp[4] = 50;
  temp[5]= 53;
  temp[6]= 59;
  temp[7] = 53;
  temp[8] = 50;
  temp[9]= 40;
  peptide_site_scores_1.push_back(temp);
  permutations.clear();
  per.clear();
  per.push_back(3);
  permutations.push_back(per);
  per.clear();
  per.push_back(6);
  permutations.push_back(per);

  ranking = ptr_test->rankWeightedPermutationPeptideScoresTest_(peptide_site_scores_1);
  ptr_test->determineHighestScoringPermutationsTest_(peptide_site_scores_1,sites,permutations,ranking);
  TEST_EQUAL(sites.size(),1)
  TEST_EQUAL(sites[0].seq_1,0)
  TEST_EQUAL(sites[0].seq_2,1)
  TEST_EQUAL(sites[0].first,3);
  TEST_EQUAL(sites[0].second,6);
  TEST_EQUAL(sites[0].peak_depth, 6)

  permutations.clear();
  per.clear();
  per.push_back(3);
  per.push_back(5);
  permutations.push_back(per);
  per.clear();
  per.push_back(5);
  per.push_back(6);
  permutations.push_back(per);
    per.clear();
  per.push_back(3);
  per.push_back(7);
  permutations.push_back(per);
    per.clear();
  per.push_back(3);
  per.push_back(6);
  permutations.push_back(per);
    per.clear();
  per.push_back(5);
  per.push_back(7);
  permutations.push_back(per);
    per.clear();
  per.push_back(6);
  per.push_back(7);
  permutations.push_back(per);
  std::vector< std::vector<double> > sec;
  peptide_site_scores_1.push_back(temp);
  peptide_site_scores_1.push_back(temp);
  peptide_site_scores_1.push_back(temp);
  peptide_site_scores_1.push_back(temp);
  
  ranking = ptr_test->rankWeightedPermutationPeptideScoresTest_(peptide_site_scores_1);
  ptr_test->determineHighestScoringPermutationsTest_(peptide_site_scores_1,sites,permutations,ranking);
  TEST_EQUAL(sites.size(),2)
  TEST_EQUAL(sites[0].seq_1,0)
  TEST_EQUAL(sites[0].seq_2,4)
  TEST_EQUAL(sites[0].first,3);
  TEST_EQUAL(sites[0].second,7);
  TEST_EQUAL(sites[0].peak_depth, 6)
  TEST_EQUAL(sites[1].seq_1,0)
  TEST_EQUAL(sites[1].seq_2,3)
  TEST_EQUAL(sites[1].first,5);
  TEST_EQUAL(sites[1].second,6);
  TEST_EQUAL(sites[1].peak_depth, 6)

}
END_SECTION

START_SECTION(computeSiteDeterminingIonsTest_(std::vector<RichPeakSpectrum> & th_spectra, ProbablePhosphoSites & candidates, std::vector<RichPeakSpectrum> & site_determining_ions))
{
  ProbablePhosphoSites candidates;
  RichPeakSpectrum temp1,temp2;
  vector<RichPeakSpectrum> site_determining_ions;
  
  RichPeakSpectrum &real_spectrum = tmp;
  std::vector<RichPeakSpectrum> windows_top10(ptr_test->peakPickingPerWindowsInSpectrumTest_(real_spectrum));
  
  AASequence seq = seq_without_phospho;
  vector<RichPeakSpectrum> th_s(ptr_test->createTheoreticalSpectraTest_(permutations, seq));
  
  candidates.seq_1 = 3;
  candidates.seq_2 = 4;
  candidates.first = 10;
  candidates.second = 7;
  ptr_test->computeSiteDeterminingIonsTest_(th_s,candidates,site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(),2)
  TEST_EQUAL(site_determining_ions[0].size(),6)
  TEST_EQUAL(site_determining_ions[1].size(),6)
  
  //=============================================================================
  
  th_s.clear();
  seq = AASequence::fromString("VTEQSP");
  candidates.seq_1 = 0;
  candidates.seq_2 = 1;
  candidates.first = 1;
  candidates.second = 4;
  
  vector<vector<Size> > p;
  perm.clear();
  perm.push_back(candidates.first);
  p.push_back(perm);

  perm.clear();
  perm.push_back(candidates.second);
  p.push_back(perm);
  
  th_s = ptr_test->createTheoreticalSpectraTest_(p, seq);
  
  ptr_test->computeSiteDeterminingIonsTest_(th_s,candidates,site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(),2)
  TEST_EQUAL(site_determining_ions[0].size(),6)
  TEST_EQUAL(site_determining_ions[1].size(),6)
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),203.102)
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),538.19)
  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),201.123)
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),540.17)
  
  candidates.first = 4;
  candidates.second = 1;
  candidates.seq_1 = 1;
  candidates.seq_2 = 0;
  
  ptr_test->computeSiteDeterminingIonsTest_(th_s,candidates,site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(),2)
  TEST_EQUAL(site_determining_ions[0].size(),6)
  TEST_EQUAL(site_determining_ions[1].size(),6)

  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),203.102)
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),538.19)
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),201.123)
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),540.17)
  
  //=============================================================================
  
  th_s.clear();
  seq = AASequence::fromString("TYQYS");
  candidates.seq_1 = 0;
  candidates.seq_2 = 1;
  candidates.first = 0;
  candidates.second = 4;
  
  p.clear();
  perm.clear();
  perm.push_back(candidates.first);
  p.push_back(perm);

  perm.clear();
  perm.push_back(candidates.second);
  p.push_back(perm);
  
  th_s = ptr_test->createTheoreticalSpectraTest_(p, seq);
  
  ptr_test->computeSiteDeterminingIonsTest_(th_s,candidates,site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(),2)
  TEST_EQUAL(site_determining_ions[0].size(),7)
  TEST_EQUAL(site_determining_ions[1].size(),7)
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),106.05)
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),636.206)
  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),186.016)
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),640.201)
  
  candidates.first = 4;
  candidates.second = 0;
  candidates.seq_1 = 1;
  candidates.seq_2 = 0;
  ptr_test->computeSiteDeterminingIonsTest_(th_s,candidates,site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(),2)
  TEST_EQUAL(site_determining_ions[0].size(),7)
  TEST_EQUAL(site_determining_ions[1].size(),7)

  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),106.05)
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),636.206)
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),186.016)
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),640.201)
  
  //=============================================================================
  
  th_s.clear();
  seq = AASequence::fromString("TSTYQYSYPP");
  candidates.seq_1 = 0;
  candidates.seq_2 = 1;
  candidates.first = 2;
  candidates.second = 6;
  
  p.clear();
  perm.clear();
  perm.push_back(candidates.first);
  p.push_back(perm);

  perm.clear();
  perm.push_back(candidates.second);
  p.push_back(perm);
  
  th_s = ptr_test->createTheoreticalSpectraTest_(p, seq);
  
  ptr_test->computeSiteDeterminingIonsTest_(th_s,candidates,site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(),2)
  TEST_EQUAL(site_determining_ions[0].size(),8)
  TEST_EQUAL(site_determining_ions[1].size(),8)

  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),370.101)
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),917.403)
  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),290.135)
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),997.37)
  
  candidates.seq_1 = 1;
  candidates.seq_2 = 0;
  candidates.first = 6;
  candidates.second = 2;
  ptr_test->computeSiteDeterminingIonsTest_(th_s,candidates,site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(),2)
  TEST_EQUAL(site_determining_ions[0].size(),8)
  TEST_EQUAL(site_determining_ions[1].size(),8)

  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),370.101)
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),917.403)
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),290.135)
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),997.37)
}
END_SECTION

START_SECTION(std::vector<Size> getSitesTest_(AASequence& without_phospho))
{
  AASequence phospho = AASequence::fromString("VTQSPSSP");
  vector<Size> tupel(ptr_test->getSitesTest_(phospho));
  TEST_EQUAL(4, tupel.size())
  TEST_EQUAL(1, tupel[0])
  TEST_EQUAL(3,tupel[1])
  TEST_EQUAL(5,tupel[2])
  TEST_EQUAL(6,tupel[3])
}
END_SECTION

START_SECTION(std::vector<std::vector<Size> > computePermutationsTest_(std::vector<Size>& tupel,Int number_of_phospho_sites))
{
  vector<Size> tupel;
  tupel.push_back(1);
  tupel.push_back(2);
  tupel.push_back(3);
  tupel.push_back(4);
  vector<vector<Size> > permutations;
  
  permutations = ptr_test->computePermutationsTest_(tupel,1);
  TEST_EQUAL(4,permutations.size())
  TEST_EQUAL(1,permutations[0][0])
  TEST_EQUAL(2,permutations[1][0])
  TEST_EQUAL(3,permutations[2][0])
  TEST_EQUAL(4,permutations[3][0])

  permutations = ptr_test->computePermutationsTest_(tupel,2);
  TEST_EQUAL(6,permutations.size())
  TEST_EQUAL(1,permutations[0][0])
  TEST_EQUAL(2,permutations[0][1])
  TEST_EQUAL(1,permutations[1][0])
  TEST_EQUAL(3,permutations[1][1])
  TEST_EQUAL(1,permutations[2][0])
  TEST_EQUAL(4,permutations[2][1])
  TEST_EQUAL(2,permutations[3][0])
  TEST_EQUAL(3,permutations[3][1])
  TEST_EQUAL(2,permutations[4][0])
  TEST_EQUAL(4,permutations[4][1])
  TEST_EQUAL(3,permutations[5][0])
  TEST_EQUAL(4,permutations[5][1])

  permutations = ptr_test->computePermutationsTest_(tupel,3);
  TEST_EQUAL(4,permutations.size())
  TEST_EQUAL(1,permutations[0][0])
  TEST_EQUAL(2,permutations[0][1])
  TEST_EQUAL(3,permutations[0][2])
  TEST_EQUAL(1,permutations[1][0])
  TEST_EQUAL(2,permutations[1][1])
  TEST_EQUAL(4,permutations[1][2])
  TEST_EQUAL(1,permutations[2][0])
  TEST_EQUAL(3,permutations[2][1])
  TEST_EQUAL(4,permutations[2][2])
  TEST_EQUAL(2,permutations[3][0])
  TEST_EQUAL(3,permutations[3][1])
  TEST_EQUAL(4,permutations[3][2])

  permutations = ptr_test->computePermutationsTest_(tupel,4);
  TEST_EQUAL(1,permutations.size())
  TEST_EQUAL(1,permutations[0][0])
  TEST_EQUAL(2,permutations[0][1])
  TEST_EQUAL(3,permutations[0][2])
  TEST_EQUAL(4,permutations[0][3])
  
  tupel.clear();
  permutations = ptr_test->computePermutationsTest_(tupel,0);
  TEST_EQUAL(0,permutations.size())
}
END_SECTION

START_SECTION(AASequence removePhosphositesFromSequenceTest_(const String sequence))
{
  String sequence = "QSSVTQVTEQS(Phospho)PK";
  TEST_EQUAL(ptr_test->removePhosphositesFromSequenceTest_(sequence).toString(),"QSSVTQVTEQSPK");
}
END_SECTION

START_SECTION(std::vector<RichPeakSpectrum> createTheoreticalSpectraTest_(const std::vector<std::vector<Size> > & permutations, const AASequence & seq_without_phospho))
{
  vector<RichPeakSpectrum> th_spectra(ptr_test->createTheoreticalSpectraTest_(permutations, seq_without_phospho));
  
  TEST_EQUAL(th_spectra.size(),5);
  TEST_EQUAL(th_spectra[0].getName(),"QS(Phospho)SVTQVTEQSPK");
  TEST_EQUAL(th_spectra[4].getName(),"QSSVTQVTEQS(Phospho)PK");
  
  th_spectra.clear();
  th_spectra = ptr_test->createTheoreticalSpectraTest_(seq_without_phospho);
  
  TEST_EQUAL(th_spectra.size(),1);
  TEST_EQUAL(th_spectra[0].getName(),seq_without_phospho.toString());
}
END_SECTION

START_SECTION(std::vector<RichPeakSpectrum> peakPickingPerWindowsInSpectrumTest_(RichPeakSpectrum & real_spectrum))
{
  // (see Beausoleil et al. Figure 3)
  MSSpectrum<RichPeak1D > tmp;  
	tmp.push_back(peak6);
  tmp.push_back(peak);
  tmp.push_back(peak2);
  tmp.push_back(peak3);
  tmp.push_back(peak4);
  tmp.push_back(peak5);
  RichPeakSpectrum &real_spectrum = tmp;
  
  std::vector<RichPeakSpectrum> windows_top10(ptr_test->peakPickingPerWindowsInSpectrumTest_(real_spectrum));
  TEST_EQUAL(windows_top10.size(),8);
  TEST_EQUAL(windows_top10[0].size(),1);
  TEST_EQUAL(windows_top10[1].size(),1);
  TEST_EQUAL(windows_top10[4].size(),0);
  TEST_EQUAL(windows_top10[7].size(),1);
}
END_SECTION

START_SECTION(Size numberOfMatchedIonsTest_(const RichPeakSpectrum & th, const RichPeakSpectrum & windows, Size depth, double fragment_mass_tolerance, bool fragment_mass_tolerance_ppm = false))
{
  RichPeakSpectrum &real_spectrum = tmp;
  double fragment_mass_tolerance = 0.5;
  bool fragment_mass_tolerance_ppm = false;
  
  vector<RichPeakSpectrum> th_spectra(ptr_test->createTheoreticalSpectraTest_(permutations, seq_without_phospho));
  std::vector<RichPeakSpectrum> windows_top10(ptr_test->peakPickingPerWindowsInSpectrumTest_(real_spectrum));
  
  //QSSVTQVTEQS(phospho)PK
  vector<RichPeakSpectrum>::iterator it = th_spectra.end() - 1;
  TEST_EQUAL(ptr_test->numberOfMatchedIonsTest_(*it, windows_top10[0], 1, fragment_mass_tolerance, fragment_mass_tolerance_ppm), 1);
  TEST_EQUAL(ptr_test->numberOfMatchedIonsTest_(*it, windows_top10[1], 1, fragment_mass_tolerance, fragment_mass_tolerance_ppm), 1);
}
END_SECTION

// of best peptide
START_SECTION(calculateCumulativeBinominalProbabilityScore)
{
  std::cout << std::endl;
  double fragment_mass_tolerance = 0.5;
  bool fragment_mass_tolerance_ppm = false;
  
  vector<ProbablePhosphoSites> phospho_sites;
  phospho_sites.clear();
  phospho_sites.resize(1);
  
  phospho_sites[0].seq_1 = 4;
  phospho_sites[0].seq_2 = 3;
  phospho_sites[0].peak_depth = 6;
  phospho_sites[0].first = 10;
  phospho_sites[0].second = 7;
  
  
  RichPeakSpectrum &real_spectrum = tmp;
  std::vector<RichPeakSpectrum> windows_top10(ptr_test->peakPickingPerWindowsInSpectrumTest_(real_spectrum));
  vector<RichPeakSpectrum> th_spectra(ptr_test->createTheoreticalSpectraTest_(permutations, seq_without_phospho));
  
  for (vector<ProbablePhosphoSites>::iterator s_it = phospho_sites.begin(); s_it < phospho_sites.end(); ++s_it)
  {
    vector<RichPeakSpectrum> site_determining_ions;
    ptr_test->computeSiteDeterminingIonsTest_(th_spectra, *s_it, site_determining_ions);
    
    Size N = site_determining_ions[0].size(); // all possibilities have the same number so take the first one
    double p = static_cast<double>(s_it->peak_depth) / 100.0;
    
    Size n_first = 0;
    for (Size depth = 0; depth != windows_top10.size(); ++depth) // for each 100 m/z window
    {
      n_first += ptr_test->numberOfMatchedIonsTest_(site_determining_ions[0], windows_top10[depth], s_it->peak_depth, fragment_mass_tolerance, fragment_mass_tolerance_ppm);
    }
    
    double P_first = ptr_test->computeCumulativeScoreTest_(N, n_first, p);
    P_first = -10 * log10(P_first);
    TEST_REAL_SIMILAR(P_first, 53.5336889240929);
  }
}
END_SECTION

delete ptr_test;
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
