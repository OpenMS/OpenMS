// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:	Petra Gutenbrunner $
// $Authors: David Wojnar, Petra Gutenbrunner $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/METADATA/DataArrays.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/AScore.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class AScoreTest : public AScore
{
  public:
    void computeSiteDeterminingIonsTest_(const std::vector<PeakSpectrum>& th_spectra, const ProbablePhosphoSites& candidates, std::vector<PeakSpectrum>& site_determining_ions) const
    {
      return computeSiteDeterminingIons_(th_spectra, candidates, site_determining_ions);
    }

    std::vector<Size> getSitesTest_(String& without_phospho) const {
      return getSites_(without_phospho);
    }

    std::vector<std::vector<Size> > computePermutationsTest_(const std::vector<Size>& sites, Int n_phosphorylation_events) const
    {
      return computePermutations_(sites, n_phosphorylation_events);
    }

    Size numberOfMatchedIonsTest_(const PeakSpectrum& th, const PeakSpectrum& windows, Size depth) const
    {
      return numberOfMatchedIons_(th, windows, depth);
    }

    //double peptideScoreTest_(const std::vector<double>& scores) const;

    void determineHighestScoringPermutationsTest_(const std::vector<std::vector<double>>& peptide_site_scores, std::vector<ProbablePhosphoSites>& sites, const std::vector<std::vector<Size>>& permutations, std::multimap<double, Size>& ranking) const
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
    
    std::vector<PeakSpectrum> createTheoreticalSpectraTest_(const std::vector<std::vector<Size>>& permutations, const AASequence& seq_without_phospho) const
    {
      return createTheoreticalSpectra_(permutations, seq_without_phospho);
    }
    
    std::vector<PeakSpectrum> peakPickingPerWindowsInSpectrumTest_(PeakSpectrum& real_spectrum) const
    {
      return peakPickingPerWindowsInSpectrum_(real_spectrum);
    }
    
    //std::vector<std::vector<double>> calculatePermutationPeptideScoresTest_(std::vector<PeakSpectrum>& th_spectra, const std::vector<PeakSpectrum>& windows_top10) const;
    
    std::multimap<double, Size> rankWeightedPermutationPeptideScoresTest_(const std::vector<std::vector<double>>& peptide_site_scores) const
    {
      return rankWeightedPermutationPeptideScores_(peptide_site_scores);
    }
};

/**
 * @brief Helper class for AScore testing with TheoreticalSpectrumGenerator
 * 
 * This class provides utility methods for generating test cases for AScore
 * using TheoreticalSpectrumGenerator to create synthetic spectra.
 */
class AScoreTestHelper
{
public:
  /**
   * @brief Generate a theoretical spectrum for a given peptide sequence
   * 
   * @param peptide The peptide sequence
   * @param intensity Intensity to set for all peaks
   * @param add_metainfo Whether to add meta information to the peaks
   * @return PeakSpectrum The generated spectrum
   */
  static PeakSpectrum generateTheoreticalSpectrum(const AASequence& peptide, 
                                                 double intensity = 100.0,
                                                 bool add_metainfo = true)
  {
    TheoreticalSpectrumGenerator spectrum_generator;
    Param params;
    params.setValue("add_metainfo", add_metainfo ? "true" : "false");
    params.setValue("add_first_prefix_ion", "true");
    params.setValue("add_losses", "false");
    spectrum_generator.setParameters(params);
    
    PeakSpectrum spectrum;
    spectrum_generator.getSpectrum(spectrum, peptide, 1, 1); // MS2 spectrum
    
    // Set all intensities to the same value for consistent testing
    for (Size i = 0; i < spectrum.size(); ++i)
    {
      spectrum[i].setIntensity(intensity);
    }
    
    // Set the name of the spectrum to the peptide sequence (important for AScore)
    spectrum.setName(peptide.toString());
    spectrum.setMSLevel(2);  // MS/MS spectrum
    spectrum.sortByPosition();
    return spectrum;
  }
  
  /**
   * @brief Generate a PeptideHit with the given sequence
   * 
   * @param sequence The peptide sequence as a string
   * @return PeptideHit The generated PeptideHit
   */
  static PeptideHit generatePeptideHit(const String& sequence)
  {
    return PeptideHit(1.0, 1, 2, AASequence::fromString(sequence));
  }
  
  /**
   * @brief Create a modified spectrum to simulate ambiguous localization
   * 
   * @param spectrum The original spectrum
   * @param ion_names Names of ions to remove (e.g., "y8", "b3")
   * @return PeakSpectrum The modified spectrum
   */
  static PeakSpectrum createAmbiguousSpectrum(const PeakSpectrum& spectrum, const std::vector<String>& ion_names)
  {
    PeakSpectrum modified_spectrum = spectrum;

    OPENMS_PRECONDITION(!modified_spectrum.getStringDataArrays().empty(), "Spectrum must have string data arrays to identify ions by name.");

    // Get the ion names data array (typically the first one)
    const auto& ion_name_array = modified_spectrum.getStringDataArrays()[0];

    // Remove peaks that would help distinguish between phosphorylation sites
    if (ion_name_array.size() == modified_spectrum.size())
    {
      std::vector<Size> indices_to_remove;
      for (Size i = 0; i < modified_spectrum.size(); ++i)
      {
        const String& ion_name = ion_name_array[i];
        for (const String& name_to_remove : ion_names)
        {
          if (ion_name.hasSubstring(name_to_remove))
          {
            indices_to_remove.push_back(i);
            break;
          }
        }
      }
      
      // Remove peaks in reverse order to avoid index shifting
      for (int i = indices_to_remove.size() - 1; i >= 0; --i)
      {
        modified_spectrum.erase(modified_spectrum.begin() + indices_to_remove[i]);
      }
      
      // Update the data arrays to match the new peak count
      for (Size i = 0; i < modified_spectrum.getStringDataArrays().size(); ++i)
      {
        auto& array = modified_spectrum.getStringDataArrays()[i];
        for (int j = indices_to_remove.size() - 1; j >= 0; --j)
        {
          array.erase(array.begin() + indices_to_remove[j]);
        }
      }
      
      for (Size i = 0; i < modified_spectrum.getIntegerDataArrays().size(); ++i)
      {
        auto& array = modified_spectrum.getIntegerDataArrays()[i];
        for (int j = indices_to_remove.size() - 1; j >= 0; --j)
        {
          array.erase(array.begin() + indices_to_remove[j]);
        }
      }
    }
    
    return modified_spectrum;
  }
};

///////////////////////////
///////////////////////////

START_TEST(AScore, "$Id$")

//=============================================================================
// create spectrum (see Beausoleil et al. Figure 3)
//=============================================================================
PeakSpectrum tmp;
DTAFile().load(OPENMS_GET_TEST_DATA_PATH("Ascore_test_input3.dta"), tmp);

//=============================================================================
  AASequence seq_without_phospho = AASequence::fromString("QSSVTQVTEQSPK");
//=============================================================================

//=============================================================================
// create permutations based on sequence QSSVTQVTEQSPK
//=============================================================================
std::vector<std::vector<Size>> permutations = { {1}, {2}, {4}, {7}, {10} };

//=============================================================================

AScore* ptr = nullptr;
AScore* nullPointer = nullptr;
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
  TEST_PRECONDITION_VIOLATED(ptr_test->computeCumulativeScoreTest_(N, n, p));

  n = 1;
  double score = ptr_test->computeCumulativeScoreTest_(N, n, p);
  TEST_REAL_SIMILAR(score, 0.1);
  N = 3;
  score = ptr_test->computeCumulativeScoreTest_(N, n, p);
  TEST_REAL_SIMILAR(score, 0.271);
}
END_SECTION

START_SECTION(determineHighestScoringPermutationsTest_(const std::vector<std::vector<double>>& peptide_site_scores, std::vector<ProbablePhosphoSites>& sites, const std::vector<std::vector<Size>>& permutations))
{
  std::multimap<double, Size> ranking;
  std::vector<std::vector<double>> peptide_site_scores_1;
  std::vector<std::vector<double>> peptide_site_scores_2;
  std::vector<std::vector<double>> peptide_site_scores_3;
  peptide_site_scores_1.resize(4);
  peptide_site_scores_2.resize(4);
  peptide_site_scores_3.resize(4);
  vector<double> temp(10, 0.1);
  peptide_site_scores_1[0] = temp;
  peptide_site_scores_2[3] = temp;
  peptide_site_scores_3[0] = temp;
  temp = vector<double>(10, 0.2);
  peptide_site_scores_1[1] = temp;
  peptide_site_scores_2[0] = temp;
  peptide_site_scores_3[3] = temp;
  temp = vector<double>(10, 0.3);
  peptide_site_scores_1[2] = temp;
  peptide_site_scores_2[1] = temp;
  peptide_site_scores_3[2] = temp;
  temp = vector<double>(10, 0.4);
  peptide_site_scores_1[3] = temp;
  peptide_site_scores_2[2] = temp;
  peptide_site_scores_3[1] = temp;

  vector<vector<Size>> permutations{ {1,3,5}, {3,5,6}, {1,3,6}, {1,5,6} };

  vector<ProbablePhosphoSites> sites;
  ranking = ptr_test->rankWeightedPermutationPeptideScoresTest_(peptide_site_scores_1);
  TEST_REAL_SIMILAR( ranking.rbegin()->first, .4);
  ptr_test->determineHighestScoringPermutationsTest_(peptide_site_scores_1, sites, permutations,ranking);
  TEST_EQUAL(sites.size(), 3);
  TEST_EQUAL(sites[0].seq_1, 3);
  TEST_EQUAL(sites[0].seq_2, 1);
  TEST_EQUAL(sites[0].second, 3);
  TEST_EQUAL(sites[0].first, 1);
  TEST_EQUAL(sites[0].peak_depth, 1);
  TEST_EQUAL(sites[1].first, 5);
  TEST_EQUAL(sites[1].second, 3);
  TEST_EQUAL(sites[1].seq_1, 3);
  TEST_EQUAL(sites[1].seq_2, 2);
  TEST_EQUAL(sites[1].peak_depth, 1);
  TEST_EQUAL(sites[2].first, 6);
  TEST_EQUAL(sites[2].second, 3);
  TEST_EQUAL(sites[2].seq_1, 3);
  TEST_EQUAL(sites[2].seq_2, 0);
  TEST_EQUAL(sites[2].peak_depth, 1);

  ranking = ptr_test->rankWeightedPermutationPeptideScoresTest_(peptide_site_scores_3);
  TEST_REAL_SIMILAR(ranking.rbegin()->first, .4);
  ptr_test->determineHighestScoringPermutationsTest_(peptide_site_scores_3, sites, permutations, ranking);
  TEST_EQUAL(sites.size(), 3);
  TEST_EQUAL(sites[0].seq_1, 1);
  TEST_EQUAL(sites[0].seq_2, 3);
  TEST_EQUAL(sites[0].second, 1);
  TEST_EQUAL(sites[0].first, 3);
  TEST_EQUAL(sites[0].peak_depth, 1);
  TEST_EQUAL(sites[1].first, 5);
  TEST_EQUAL(sites[1].second, 1);
  TEST_EQUAL(sites[1].seq_1, 1);
  TEST_EQUAL(sites[1].seq_2, 2);
  TEST_EQUAL(sites[1].peak_depth, 1);
  TEST_EQUAL(sites[2].first, 6);
  TEST_EQUAL(sites[2].second, 1);
  TEST_EQUAL(sites[2].seq_1, 1);
  TEST_EQUAL(sites[2].seq_2, 0);
  TEST_EQUAL(sites[2].peak_depth, 1);

  ranking = ptr_test->rankWeightedPermutationPeptideScoresTest_(peptide_site_scores_2);
  TEST_REAL_SIMILAR(ranking.rbegin()->first, .4);
  ptr_test->determineHighestScoringPermutationsTest_(peptide_site_scores_2, sites, permutations, ranking);
  TEST_EQUAL(sites.size(), 3);
  TEST_EQUAL(sites[0].seq_1, 2);
  TEST_EQUAL(sites[0].seq_2, 1);
  TEST_EQUAL(sites[0].second, 5);
  TEST_EQUAL(sites[0].first, 1);
  TEST_EQUAL(sites[0].peak_depth, 1);
  TEST_EQUAL(sites[1].first, 3);
  TEST_EQUAL(sites[1].second, 5);
  TEST_EQUAL(sites[1].seq_1, 2);
  TEST_EQUAL(sites[1].seq_2, 3);
  TEST_EQUAL(sites[1].peak_depth, 1);
  TEST_EQUAL(sites[2].first, 6);
  TEST_EQUAL(sites[2].second, 5);
  TEST_EQUAL(sites[2].seq_1, 2);
  TEST_EQUAL(sites[2].seq_2, 0);
  TEST_EQUAL(sites[2].peak_depth, 1)

  peptide_site_scores_1.clear();
  temp = {55, 60, 75, 100, 90, 120, 125, 120, 100, 90};
  peptide_site_scores_1.push_back(temp);
  temp = { 40, 50, 53, 60, 50, 53, 59, 53, 50, 40 };
  peptide_site_scores_1.push_back(temp);
  permutations = { {3}, {6} };

  ranking = ptr_test->rankWeightedPermutationPeptideScoresTest_(peptide_site_scores_1);
  TEST_REAL_SIMILAR(ranking.rbegin()->first, 94.10714285714287);
  ptr_test->determineHighestScoringPermutationsTest_(peptide_site_scores_1, sites, permutations, ranking);
  TEST_EQUAL(sites.size(), 1)
  TEST_EQUAL(sites[0].seq_1, 0)
  TEST_EQUAL(sites[0].seq_2, 1)
  TEST_EQUAL(sites[0].first, 3);
  TEST_EQUAL(sites[0].second, 6);
  TEST_EQUAL(sites[0].peak_depth, 6)

  permutations = { {3, 5}, {5, 6}, {3, 7}, {3, 6}, {5, 7}, {6, 7} };

  peptide_site_scores_1.push_back(temp);
  peptide_site_scores_1.push_back(temp);
  peptide_site_scores_1.push_back(temp);
  peptide_site_scores_1.push_back(temp);
  
  ranking = ptr_test->rankWeightedPermutationPeptideScoresTest_(peptide_site_scores_1);
  ptr_test->determineHighestScoringPermutationsTest_(peptide_site_scores_1, sites, permutations, ranking);
  TEST_EQUAL(sites.size(), 2);
  TEST_EQUAL(sites[0].seq_1, 0);
  TEST_EQUAL(sites[0].seq_2, 4);
  TEST_EQUAL(sites[0].first, 3);
  TEST_EQUAL(sites[0].second, 7);
  TEST_EQUAL(sites[0].peak_depth, 6);
  TEST_EQUAL(sites[1].seq_1, 0);
  TEST_EQUAL(sites[1].seq_2, 3);
  TEST_EQUAL(sites[1].first, 5);
  TEST_EQUAL(sites[1].second, 6);
  TEST_EQUAL(sites[1].peak_depth, 6);
}
END_SECTION

START_SECTION(computeSiteDeterminingIonsTest_(const std::vector<PeakSpectrum>& th_spectra, const ProbablePhosphoSites& candidates, std::vector<PeakSpectrum>& site_determining_ions))
{
  ProbablePhosphoSites candidates;
  PeakSpectrum temp1, temp2;
  vector<PeakSpectrum> site_determining_ions;
  
  AASequence seq = seq_without_phospho;
  vector<PeakSpectrum> th_s = ptr_test->createTheoreticalSpectraTest_(permutations, seq);
  
  candidates.seq_1 = 3;
  candidates.seq_2 = 4;
  candidates.first = 10;
  candidates.second = 7;
  ptr_test->computeSiteDeterminingIonsTest_(th_s, candidates, site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(), 2);
  TEST_EQUAL(site_determining_ions[0].size(), 6);
  TEST_EQUAL(site_determining_ions[1].size(), 6);
  
  //=============================================================================
  
  th_s.clear();
  seq = AASequence::fromString("VTEQSP");
  candidates.seq_1 = 0;
  candidates.seq_2 = 1;
  candidates.first = 1;
  candidates.second = 4;
  
  vector<vector<Size>> p { {candidates.first}, {candidates.second} };
  
  th_s = ptr_test->createTheoreticalSpectraTest_(p, seq);
  
  ptr_test->computeSiteDeterminingIonsTest_(th_s, candidates, site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(), 2);
  TEST_EQUAL(site_determining_ions[0].size(), 6);
  TEST_EQUAL(site_determining_ions[1].size(), 6);
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(), 203.102);
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size() - 1].getMZ(), 538.19);
  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(), 201.123);
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size() - 1].getMZ(), 540.17);
  
  candidates.first = 4;
  candidates.second = 1;
  candidates.seq_1 = 1;
  candidates.seq_2 = 0;
  
  ptr_test->computeSiteDeterminingIonsTest_(th_s, candidates, site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(), 2);
  TEST_EQUAL(site_determining_ions[0].size(), 6);
  TEST_EQUAL(site_determining_ions[1].size(), 6);

  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(), 203.102);
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size() - 1].getMZ(), 538.19);
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(), 201.123);
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size() - 1].getMZ(), 540.17);
  
  //=============================================================================
  
  th_s.clear();
  seq = AASequence::fromString("TYQYS");
  candidates.seq_1 = 0;
  candidates.seq_2 = 1;
  candidates.first = 0;
  candidates.second = 4;
  
  p = { { candidates.first },{ candidates.second } };

  
  th_s = ptr_test->createTheoreticalSpectraTest_(p, seq);
  
  ptr_test->computeSiteDeterminingIonsTest_(th_s, candidates, site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(), 2);
  TEST_EQUAL(site_determining_ions[0].size(), 8);
  TEST_EQUAL(site_determining_ions[1].size(), 8);
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(), 106.05);
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size() - 1].getMZ(), 636.206);
  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(), 102.0550);
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size() - 1].getMZ(), 640.201);
  
  candidates.first = 4;
  candidates.second = 0;
  candidates.seq_1 = 1;
  candidates.seq_2 = 0;
  ptr_test->computeSiteDeterminingIonsTest_(th_s, candidates, site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(), 2);
  TEST_EQUAL(site_determining_ions[0].size(), 8);
  TEST_EQUAL(site_determining_ions[1].size(), 8);

  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(), 106.05);
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size() - 1].getMZ(), 636.206);
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(), 102.0550);
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size() - 1].getMZ(), 640.201);
  
  //=============================================================================
  
  th_s.clear();
  seq = AASequence::fromString("TSTYQYSYPP");
  candidates.seq_1 = 0;
  candidates.seq_2 = 1;
  candidates.first = 2;
  candidates.second = 6;

  p = { { candidates.first },{ candidates.second } };
  
  th_s = ptr_test->createTheoreticalSpectraTest_(p, seq);
  
  ptr_test->computeSiteDeterminingIonsTest_(th_s, candidates, site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(), 2);
  TEST_EQUAL(site_determining_ions[0].size(), 8);
  TEST_EQUAL(site_determining_ions[1].size(), 8);

  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(), 370.101);
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size() - 1].getMZ(), 917.403);
  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(), 290.135);
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size() - 1].getMZ(), 997.37);
  
  candidates.seq_1 = 1;
  candidates.seq_2 = 0;
  candidates.first = 6;
  candidates.second = 2;
  ptr_test->computeSiteDeterminingIonsTest_(th_s, candidates, site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(), 2);
  TEST_EQUAL(site_determining_ions[0].size(), 8);
  TEST_EQUAL(site_determining_ions[1].size(), 8);

  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(), 370.101);
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size() - 1].getMZ(), 917.403);
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(), 290.135);
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size() - 1].getMZ(), 997.37);
  
  //=============================================================================
  
  //ATPGNLGSSVLMY(Phospho)K; ATPGNLGSS(Phospho)VLMYK
  th_s.clear();
  seq = AASequence::fromString("ATPGNLGSSVLMYK");
  candidates.seq_1 = 0;
  candidates.seq_2 = 1;
  candidates.first = 12;
  candidates.second = 8;
  
  p = { { candidates.first },{ candidates.second } };
  
  th_s = ptr_test->createTheoreticalSpectraTest_(p, seq);
  
  ptr_test->computeSiteDeterminingIonsTest_(th_s, candidates, site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(), 2);
  TEST_EQUAL(site_determining_ions[0].size(), 8);
  TEST_EQUAL(site_determining_ions[1].size(), 4);
  
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(), 390.142);
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size() - 1].getMZ(), 1128.57);
  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(), 310.176);
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size() - 1].getMZ(), 1208.54);
    
  candidates.seq_1 = 1;
  candidates.seq_2 = 0;
  candidates.first = 8;
  candidates.second = 12;
  ptr_test->computeSiteDeterminingIonsTest_(th_s, candidates, site_determining_ions);
  TEST_EQUAL(site_determining_ions.size(), 2);
  TEST_EQUAL(site_determining_ions[0].size(), 4);
  TEST_EQUAL(site_determining_ions[1].size(), 8);

  TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(), 390.142);
  TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size() - 1].getMZ(), 1128.57);
  TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(), 310.176);
  TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size() - 1].getMZ(), 1208.54);
}
END_SECTION

START_SECTION(bool isPhosphoDecoy(char residue))
{
  // Test PhosphoDecoy sites recognition
  TEST_EQUAL(ptr_test->isPhosphoDecoySite('A'), true)
  TEST_EQUAL(ptr_test->isPhosphoDecoySite('G'), true)
  TEST_EQUAL(ptr_test->isPhosphoDecoySite('L'), true)
  
  // Test non-PhosphoDecoy sites
  TEST_EQUAL(ptr_test->isPhosphoDecoySite('S'), false)
  TEST_EQUAL(ptr_test->isPhosphoDecoySite('T'), false)
  TEST_EQUAL(ptr_test->isPhosphoDecoySite('Y'), false)
}
END_SECTION

START_SECTION(std::vector<Size> getSitesTest_(const AASequence& without_phospho))
{
  AASequence phospho = AASequence::fromString("VTQSPSSP");
  String unmodified_sequence = phospho.toUnmodifiedString();
  vector<Size> tupel(ptr_test->getSitesTest_(unmodified_sequence));
  TEST_EQUAL(4, tupel.size()) // vTqSpSSp
  TEST_EQUAL(1, tupel[0])
  TEST_EQUAL(3, tupel[1])
  TEST_EQUAL(5, tupel[2])
  TEST_EQUAL(6, tupel[3])
  
  // Test with add_decoys=true
  Param params;
  params.setValue("add_decoys", "true");
  ptr_test->setParameters(params);
  
  AASequence phospho_decoy = AASequence::fromString("VTLAGSPS");
  String unmodified_sequence_decoy = phospho_decoy.toUnmodifiedString();
  vector<Size> tupel_with_decoy(ptr_test->getSitesTest_(unmodified_sequence_decoy));
  TEST_EQUAL(6, tupel_with_decoy.size())  // vTLAGSpS
  TEST_EQUAL(tupel_with_decoy[0], 1) // T
  TEST_EQUAL(tupel_with_decoy[1], 2) // L
  TEST_EQUAL(tupel_with_decoy[2], 3) // A
  TEST_EQUAL(tupel_with_decoy[3], 4) // G
  TEST_EQUAL(tupel_with_decoy[4], 5) // S
  TEST_EQUAL(tupel_with_decoy[5], 7) // S
  
  // Test with add_decoys=false
  params.setValue("add_decoys", "false");
  ptr_test->setParameters(params);
  vector<Size> tupel_no_decoy(ptr_test->getSitesTest_(unmodified_sequence_decoy));
  TEST_EQUAL(tupel_no_decoy.size(), 3) // S and P residues should be found
}
END_SECTION

START_SECTION(std::vector<std::vector<Size>> computePermutationsTest_(const std::vector<Size>& tupel, Int number_of_phospho_sites))
{
  vector<Size> tupel{1, 2, 3, 4};
  vector<vector<Size> > permutations;
  
  permutations = ptr_test->computePermutationsTest_(tupel, 1);
  TEST_EQUAL(4, permutations.size());
  TEST_EQUAL(1, permutations[0][0]);
  TEST_EQUAL(2, permutations[1][0]);
  TEST_EQUAL(3, permutations[2][0]);
  TEST_EQUAL(4, permutations[3][0]);

  permutations = ptr_test->computePermutationsTest_(tupel, 2);
  TEST_EQUAL(6, permutations.size());
  TEST_EQUAL(1, permutations[0][0]);
  TEST_EQUAL(2, permutations[0][1]);
  TEST_EQUAL(1, permutations[1][0]);
  TEST_EQUAL(3, permutations[1][1]);
  TEST_EQUAL(1, permutations[2][0]);
  TEST_EQUAL(4, permutations[2][1]);
  TEST_EQUAL(2, permutations[3][0]);
  TEST_EQUAL(3, permutations[3][1]);
  TEST_EQUAL(2, permutations[4][0]);
  TEST_EQUAL(4, permutations[4][1]);
  TEST_EQUAL(3, permutations[5][0]);
  TEST_EQUAL(4, permutations[5][1]);

  permutations = ptr_test->computePermutationsTest_(tupel, 3);
  TEST_EQUAL(4, permutations.size());
  TEST_EQUAL(1, permutations[0][0]);
  TEST_EQUAL(2, permutations[0][1]);
  TEST_EQUAL(3, permutations[0][2]);
  TEST_EQUAL(1, permutations[1][0]);
  TEST_EQUAL(2, permutations[1][1]);
  TEST_EQUAL(4, permutations[1][2]);
  TEST_EQUAL(1, permutations[2][0]);
  TEST_EQUAL(3, permutations[2][1]);
  TEST_EQUAL(4, permutations[2][2]);
  TEST_EQUAL(2, permutations[3][0]);
  TEST_EQUAL(3, permutations[3][1]);
  TEST_EQUAL(4, permutations[3][2]);

  permutations = ptr_test->computePermutationsTest_(tupel, 4);
  TEST_EQUAL(1, permutations.size());
  TEST_EQUAL(1, permutations[0][0]);
  TEST_EQUAL(2, permutations[0][1]);
  TEST_EQUAL(3, permutations[0][2]);
  TEST_EQUAL(4, permutations[0][3]);
  
  tupel.clear();
  permutations = ptr_test->computePermutationsTest_(tupel, 0);
  TEST_EQUAL(0, permutations.size());
}
END_SECTION

START_SECTION(AASequence removePhosphositesFromSequenceTest_(const String sequence))
{
  String sequence = "QSSVTQVTEQS(Phospho)PK";
  String sequence_decoy = "QA(PhosphoDecoy)SVTQVTEQS(Phospho)PK";
  
  TEST_EQUAL(ptr_test->removePhosphositesFromSequenceTest_(sequence).toString(), "QSSVTQVTEQSPK");
  TEST_EQUAL(ptr_test->removePhosphositesFromSequenceTest_(sequence_decoy).toString(), "QASVTQVTEQSPK");
}
END_SECTION

START_SECTION(std::vector<PeakSpectrum> peakPickingPerWindowsInSpectrumTest_(PeakSpectrum& real_spectrum))
{
  PeakSpectrum& real_spectrum = tmp;
  
  std::vector<PeakSpectrum> windows_top10 = ptr_test->peakPickingPerWindowsInSpectrumTest_(real_spectrum);
  TEST_EQUAL(windows_top10.size(), 8);
  TEST_EQUAL(windows_top10[0].size(), 1);
  TEST_EQUAL(windows_top10[1].size(), 1);
  TEST_EQUAL(windows_top10[4].size(), 0);
  TEST_EQUAL(windows_top10[7].size(), 1);
}
END_SECTION

START_SECTION(Size numberOfMatchedIonsTest_(const PeakSpectrum& th, const PeakSpectrum& windows, Size depth))
{
  PeakSpectrum& real_spectrum = tmp;
  Param params;
  params.setValue("fragment_mass_tolerance", 0.5);
  ptr_test->setParameters(params);
 
  vector<PeakSpectrum> th_spectra = ptr_test->createTheoreticalSpectraTest_(permutations, seq_without_phospho);
  std::vector<PeakSpectrum> windows_top10 = ptr_test->peakPickingPerWindowsInSpectrumTest_(real_spectrum);
  
  //QSSVTQVTEQS(phospho)PK
  vector<PeakSpectrum>::iterator it = th_spectra.end() - 1;
  TEST_EQUAL(ptr_test->numberOfMatchedIonsTest_(*it, windows_top10[0], 1), 1);
  TEST_EQUAL(ptr_test->numberOfMatchedIonsTest_(*it, windows_top10[1], 1), 1);
}
END_SECTION

// of best peptide
START_SECTION(calculateCumulativeBinominalProbabilityScore)
{
  vector<ProbablePhosphoSites> phospho_sites;
  phospho_sites.clear();
  phospho_sites.resize(1);
  
  phospho_sites[0].seq_1 = 4;
  phospho_sites[0].seq_2 = 3;
  phospho_sites[0].peak_depth = 6;
  phospho_sites[0].first = 10;
  phospho_sites[0].second = 7;
  
  
  PeakSpectrum& real_spectrum = tmp;
  std::vector<PeakSpectrum> windows_top10 = ptr_test->peakPickingPerWindowsInSpectrumTest_(real_spectrum);
  vector<PeakSpectrum> th_spectra = ptr_test->createTheoreticalSpectraTest_(permutations, seq_without_phospho);
  
  for (vector<ProbablePhosphoSites>::iterator s_it = phospho_sites.begin(); s_it < phospho_sites.end(); ++s_it)
  {
    vector<PeakSpectrum> site_determining_ions;
    ptr_test->computeSiteDeterminingIonsTest_(th_spectra, *s_it, site_determining_ions);
    
    Size N = site_determining_ions[0].size(); // all possibilities have the same number so take the first one
    double p = static_cast<double>(s_it->peak_depth) / 100.0;
    
    Size n_first = 0;
    for (Size depth = 0; depth != windows_top10.size(); ++depth) // for each 100 m/z window
    {
      n_first += ptr_test->numberOfMatchedIonsTest_(site_determining_ions[0], windows_top10[depth], s_it->peak_depth);
    }

    double P_first = ptr_test->computeCumulativeScoreTest_(N, n_first, p);
    P_first = -10 * log10(P_first);
    TEST_REAL_SIMILAR(P_first, 53.5336889240929);
  }
}
END_SECTION

START_SECTION(std::vector<PeakSpectrum> createTheoreticalSpectraTest_(const std::vector<std::vector<Size>>& permutations, const AASequence& seq_without_phospho))
{
  // create theoretical based on permutations
  vector<PeakSpectrum> th_spectra(ptr_test->createTheoreticalSpectraTest_(permutations, seq_without_phospho));
  TEST_EQUAL(th_spectra.size(), 5);
  TEST_EQUAL(th_spectra[0].getName(), "QS(Phospho)SVTQVTEQSPK");
  TEST_EQUAL(th_spectra[4].getName(), "QSSVTQVTEQS(Phospho)PK");
  TEST_REAL_SIMILAR(th_spectra[4][0].getMZ(), 129.0658);
  TEST_REAL_SIMILAR(th_spectra[4][2].getMZ(), 216.0979);
  TEST_REAL_SIMILAR(th_spectra[4][21].getMZ(), 1283.5879);
  
  th_spectra.clear();
}
END_SECTION 

START_SECTION(PeptideHit AScore::compute(const PeptideHit& hit, PeakSpectrum& real_spectrum) const)
{
  // ====================================================================================================================================
  // The Ascore results differ to the results of the Ascore tool provided on the website http://ascore.med.harvard.edu/ascore.html
  // But it seems that the online version has some issues calculating the Ascore using the cumulative binomial probability formula.
  // E.g. with the values 6, 5, 0.06 for the variables N, n, p the calculated Ascore using WolframAlpha is 53.5337, which does not 
  // conform to the result 53.57, which is mentioned in the paper (see Fig. 3c).
  // In addition the site determining ions calculation seems not reliable, because in some test cases more site determining ions 
  // were calculated than it could be possible.
  // Another reason for the differences of the results could be the fragment ion tolerance used to match the theoretical spectra
  // with the real spectra. The value used in the Ascore tool provided on the website is not mentioned.  
  // ====================================================================================================================================
  
  PeakSpectrum real_spectrum;
  Param params;
  params.setValue("fragment_mass_tolerance", 0.6);
  params.setValue("add_decoys", "false");
  ptr_test->setParameters(params);
  
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("Ascore_test_input1.dta"), real_spectrum);
  PeptideHit hit1(1.0, 1, 1, AASequence::fromString("QSSVT(Phospho)QSK"));
  hit1 = ptr_test->compute(hit1, real_spectrum);
  
  // http://ascore.med.harvard.edu/ascore.html result=3.51, sequence=QSSVT*QSK
  TEST_REAL_SIMILAR(static_cast<double>(hit1.getMetaValue("AScore_1")), 8.65157151899052);
  TEST_EQUAL(hit1.getSequence().toString(), "QSS(Phospho)VTQSK");
  
  // ===========================================================================
  
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("Ascore_test_input2.dta"), real_spectrum);
  PeptideHit hit2(1.0, 1, 1, AASequence::fromString("RIRLT(Phospho)ATTR"));
  hit2 = ptr_test->compute(hit2, real_spectrum);
  
  // http://ascore.med.harvard.edu/ascore.html result=21.3
  TEST_REAL_SIMILAR(static_cast<double>(hit2.getMetaValue("AScore_1")), 18.8755623850511);
  TEST_EQUAL(hit2.getSequence().toString(), "RIRLT(Phospho)ATTR");
  
  // ===========================================================================
  
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("Ascore_test_input3.dta"), real_spectrum);
  PeptideHit hit3(1.0, 1, 1, AASequence::fromString("QSSVTQVTEQS(Phospho)PK"));
  hit3 = ptr_test->compute(hit3, real_spectrum);
  
  // http://ascore.med.harvard.edu/ascore.html result=88.3
  TEST_REAL_SIMILAR(static_cast<double>(hit3.getMetaValue("AScore_1")), 88.3030731386678);
  TEST_EQUAL(hit3.getSequence().toString(), "QSSVTQVTEQS(Phospho)PK"); 
  
  // ===========================================================================
  
  params.setValue("fragment_mass_tolerance", 0.05);
  ptr_test->setParameters(params);
  
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("Ascore_test_input4.dta"), real_spectrum);
  PeptideHit hit4(1.0, 1, 1, AASequence::fromString("ATPGNLGSSVLHS(Phospho)K"));
  
  hit4 = ptr_test->compute(hit4, real_spectrum);
  
  // http://ascore.med.harvard.edu/ascore.html result=88.3
  TEST_REAL_SIMILAR(static_cast<double>(hit4.getMetaValue("AScore_1")), 49.2714597801023);
  TEST_EQUAL(hit4.getSequence().toString(), "ATPGNLGSSVLHS(Phospho)K");
  
  // ===========================================================================
  // PPM UNIT TEST
  // ===========================================================================
  
  params.setValue("fragment_mass_tolerance", 700.0); // 0.6 Da were converted to ppm based on a small peptide 
  params.setValue("fragment_mass_unit", "ppm");
  
  ptr_test->setParameters(params);
  
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("Ascore_test_input1.dta"), real_spectrum);
  PeptideHit hit5(1.0, 1, 1, AASequence::fromString("QSSVT(Phospho)QSK"));
  hit5 = ptr_test->compute(hit5, real_spectrum);
  
  // http://ascore.med.harvard.edu/ascore.html result=3.51, sequence=QSSVT*QSK
  TEST_REAL_SIMILAR(static_cast<double>(hit5.getMetaValue("AScore_1")), 6.53833235677545);
  TEST_EQUAL(hit5.getSequence().toString(), "QSS(Phospho)VTQSK");
  
  params.setValue("fragment_mass_tolerance", 70.0); // 0.05 Da were converted to ppm based on a small peptide
  ptr_test->setParameters(params);
  
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("Ascore_test_input4.dta"), real_spectrum);
  PeptideHit hit6(1.0, 1, 1, AASequence::fromString("ATPGNLGSSVLHS(Phospho)K"));
  
  hit6 = ptr_test->compute(hit6, real_spectrum);
  
  // http://ascore.med.harvard.edu/ascore.html result=88.3
  TEST_REAL_SIMILAR(static_cast<double>(hit6.getMetaValue("AScore_1")), 40.6506162613816);
  TEST_EQUAL(hit6.getSequence().toString(), "ATPGNLGSSVLHS(Phospho)K");

  // ===========================================================================
  // check if special score is used for unambiguous assignment:
  PeptideHit hit7(1.0, 1, 1, AASequence::fromString("PEPT(Phospho)IDE"));
  hit7 = ptr_test->compute(hit7, real_spectrum);
  TEST_REAL_SIMILAR(hit7.getScore(), 1e6);
  
  // ===========================================================================
  // Test PhosphoDecoy functionality
  // ===========================================================================
  
  params.setValue("add_decoys", "true");
  ptr_test->setParameters(params);
  
  // Test with PhosphoDecoy modification
  PeptideHit hit_decoy(1.0, 1, 1, AASequence::fromString("ATPL(PhosphoDecoy)NLGSSVLHSK"));
  hit_decoy = ptr_test->compute(hit_decoy, real_spectrum);
  
  // Check that metadata is stored correctly
  TEST_EQUAL(hit_decoy.getMetaValue("regular_phospho_count"), 0);
  TEST_EQUAL(hit_decoy.getMetaValue("phospho_decoy_count"), 1);
  
  // Test with both regular phosphorylation and PhosphoDecoy
  PeptideHit hit_mixed(1.0, 1, 1, AASequence::fromString("ATPL(PhosphoDecoy)NLGSSVLHS(Phospho)K"));
  hit_mixed = ptr_test->compute(hit_mixed, real_spectrum);
  
  // Check that both types are recognized and counted
  TEST_EQUAL(hit_mixed.getMetaValue("regular_phospho_count"), 1);
  TEST_EQUAL(hit_mixed.getMetaValue("phospho_decoy_count"), 1);
  
  // Test with PhosphoDecoy disabled
  params.setValue("add_decoys", "false");
  ptr_test->setParameters(params);
  
  PeptideHit hit_disabled(1.0, 1, 1, AASequence::fromString("ATPL(PhosphoDecoy)NLGSSVLHSK"));
  hit_disabled = ptr_test->compute(hit_disabled, real_spectrum);
  
  // The PhosphoDecoy sites should be ignored when the option is disabled
  TEST_EQUAL(hit_disabled.metaValueExists("phospho_decoy_count"), false);
  
  // ===========================================================================
  // Test ProForma output
  // ===========================================================================
  
  // Test with regular phosphorylation
  params.setValue("add_decoys", "false");
  ptr_test->setParameters(params);
  
  PeptideHit hit_proforma(1.0, 1, 1, AASequence::fromString("AT(Phospho)PGNLGSSVLHS(Phospho)K"));
  hit_proforma = ptr_test->compute(hit_proforma, real_spectrum);
  
  // Check that ProForma meta value exists
  TEST_EQUAL(hit_proforma.metaValueExists("ProForma"), true);
  
  // Check the format of the ProForma
  String proforma_str = hit_proforma.getMetaValue("ProForma");
  TEST_EQUAL(proforma_str, "AT[Phospho|score=0.9567]PGNLGSSVLHS[Phospho|score=0.9760]K");
  
  // Test with both regular phosphorylation and PhosphoDecoy
  params.setValue("add_decoys", "true");
  ptr_test->setParameters(params);
  
  PeptideHit hit_proforma_mixed(1.0, 1, 1, AASequence::fromString("ATPL(PhosphoDecoy)NLGSSVLHS(Phospho)K"));
  hit_proforma_mixed = ptr_test->compute(hit_proforma_mixed, real_spectrum);
  
  // Check that ProForma meta value exists and contains both modification types
  TEST_EQUAL(hit_proforma_mixed.metaValueExists("ProForma"), true);
  proforma_str = hit_proforma_mixed.getMetaValue("ProForma");
  TEST_EQUAL(proforma_str, "ATPL[PhosphoDecoy|score=0.9239]NLGSSVLHS[Phospho|score=0.9718]K");
}
END_SECTION 

START_SECTION(Enhanced AScore Tests Using TheoreticalSpectrumGenerator)
{
  // Create AScore instance with default parameters
  AScore ascore;
  Param params;
  params.setValue("fragment_mass_tolerance", 0.1);
  params.setValue("fragment_mass_unit", "Da");
  params.setValue("add_decoys", "false");
  ascore.setParameters(params);
  
  // Test case 1: Single phosphorylation site with clear localization
  {
    // Create peptide with one phosphorylation site
    AASequence peptide = AASequence::fromString("PEPTIDES(Phospho)");
    
    // Generate theoretical spectrum
    PeakSpectrum spectrum = AScoreTestHelper::generateTheoreticalSpectrum(peptide);
    
    // Create PeptideHit with the same sequence
    PeptideHit hit = AScoreTestHelper::generatePeptideHit("PEPTIDES(Phospho)");
    
    // Run AScore
    hit = ascore.compute(hit, spectrum);
    
    // Verify results
    TEST_EQUAL(hit.getSequence().toString(), "PEPTIDES(Phospho)");
    TEST_EQUAL(hit.metaValueExists("ProForma"), true);
    
    // Verify ProForma string
    String proforma = hit.getMetaValue("ProForma");
    TEST_EQUAL(proforma, "PEPTIDES[Phospho|score=0.9948]");
  }
  
  // Test case 2: Multiple phosphorylation sites
  {
    // Create peptide with multiple phosphorylation sites
    AASequence peptide = AASequence::fromString("PEPT(Phospho)IDES(Phospho)");
    PeptideHit hit = AScoreTestHelper::generatePeptideHit(peptide.toString());
    TEST_EQUAL(hit.getSequence().toString(), "PEPT(Phospho)IDES(Phospho)");

    // Generate theoretical spectrum
    PeakSpectrum spectrum = AScoreTestHelper::generateTheoreticalSpectrum(peptide);
    TEST_TRUE(spectrum.size() > 0);

    // Run AScore
    hit = ascore.compute(hit, spectrum);
    
    // Verify results
    TEST_TRUE(hit.metaValueExists("ProForma"));
    
    // Verify ProForma string contains both sites
    String proforma = hit.getMetaValue("ProForma");
    TEST_EQUAL(proforma, "PEPT[Phospho|score=1.0000]IDES[Phospho|score=1.0000]");
  }
  
  // Test case 3: PhosphoDecoy sites
  {
    // Enable PhosphoDecoy
    params.setValue("add_decoys", "true");
    ascore.setParameters(params);
    
    // Create peptide with PhosphoDecoy site
    AASequence peptide = AASequence::fromString("PEPTIDEA(PhosphoDecoy)");
    
    // Generate theoretical spectrum
    PeakSpectrum spectrum = AScoreTestHelper::generateTheoreticalSpectrum(peptide);
    
    // Create PeptideHit with the same sequence
    PeptideHit hit = AScoreTestHelper::generatePeptideHit("PEPTIDEA(PhosphoDecoy)");
    
    // Run AScore
    hit = ascore.compute(hit, spectrum);
    
    // Verify results
    TEST_EQUAL(hit.getSequence().toString(), "PEPTIDEA(PhosphoDecoy)");
    TEST_TRUE(hit.metaValueExists("ProForma"));
    
    // Verify ProForma string
    String proforma = hit.getMetaValue("ProForma");
    TEST_EQUAL(proforma, "PEPTIDEA[PhosphoDecoy|score=0.9948]");
    
    // Verify metadata
    TEST_EQUAL(hit.getMetaValue("phospho_decoy_count"), 1);
    TEST_EQUAL(hit.getMetaValue("regular_phospho_count"), 0);
  }
  
  // Test case 4: Mix of Phospho and PhosphoDecoy sites
  {
    // Create peptide with both types of sites
    AASequence peptide = AASequence::fromString("PEPT(Phospho)IDEA(PhosphoDecoy)");
    
    // Generate theoretical spectrum
    PeakSpectrum spectrum = AScoreTestHelper::generateTheoreticalSpectrum(peptide);
    
    // Create PeptideHit with the same sequence
    PeptideHit hit = AScoreTestHelper::generatePeptideHit("PEPT(Phospho)IDEA(PhosphoDecoy)");
    
    // Run AScore
    hit = ascore.compute(hit, spectrum);
    
    // Verify results
    TEST_EQUAL(hit.getSequence().toString(), "PEPT(Phospho)IDEA(PhosphoDecoy)");
    TEST_EQUAL(hit.metaValueExists("ProForma"), true);
    
    // Verify ProForma string contains both modification types
    String proforma = hit.getMetaValue("ProForma");
    TEST_EQUAL(proforma, "PEPT[Phospho|score=1.0000]IDEA[PhosphoDecoy|score=1.0000]");
    
    // Verify metadata
    TEST_EQUAL(hit.getMetaValue("phospho_decoy_count"), 1);
    TEST_EQUAL(hit.getMetaValue("regular_phospho_count"), 1);
  }
  
  // Test case 5: All sites phosphorylated (unambiguous assignment)
  {
    // Create peptide with all S/T/Y residues phosphorylated
    AASequence peptide = AASequence::fromString("PEPT(Phospho)IDES(Phospho)");
    
    // Generate theoretical spectrum
    PeakSpectrum spectrum = AScoreTestHelper::generateTheoreticalSpectrum(peptide);
    
    // Create PeptideHit with the same sequence
    PeptideHit hit = AScoreTestHelper::generatePeptideHit("PEPT(Phospho)IDES(Phospho)");
    
    // Run AScore
    hit = ascore.compute(hit, spectrum);
    
    // Verify results - expect unambiguous score and ProForma string contains all sites
    String proforma = hit.getMetaValue("ProForma");
    TEST_EQUAL(proforma, "PEPT[Phospho|score=1.0000]IDES[Phospho|score=1.0000]");
  }
  
  // Test case 6: Ambiguous phosphorylation site localization
  {
    // Create a peptide with multiple potential phosphorylation sites
    // but only one actual phosphorylation
    AASequence peptide = AASequence::fromString("PEPTIDES(Phospho)T");
    
    // Generate theoretical spectrum
    PeakSpectrum spectrum = AScoreTestHelper::generateTheoreticalSpectrum(peptide);
    
    // Create a modified spectrum with some site-determining ions removed
    // to make localization ambiguous
    std::vector<String> ions_to_remove = {"y2", "y3", "b7", "b8"};
    PeakSpectrum ambiguous_spectrum = AScoreTestHelper::createAmbiguousSpectrum(spectrum, ions_to_remove);
    
    // Create PeptideHit with the same sequence
    PeptideHit hit = AScoreTestHelper::generatePeptideHit("PEPTIDES(Phospho)T");
    
    // Run AScore with the ambiguous spectrum
    hit = ascore.compute(hit, ambiguous_spectrum);
    
    // Verify ProForma string exists
    TEST_EQUAL(hit.metaValueExists("ProForma"), true);
    
    // The localization might be ambiguous, but we can't predict the exact outcome
    // since it depends on the specific ions that were removed
  }
  
  // Test case 7: Peptide with many potential phosphorylation sites but only a few actual phosphorylations
  {
    // Create a peptide with many S/T/Y residues but only some phosphorylated
    AASequence peptide = AASequence::fromString("STYSTSYSTSTYS(Phospho)TSYT(Phospho)");
    
    // Generate theoretical spectrum
    PeakSpectrum spectrum = AScoreTestHelper::generateTheoreticalSpectrum(peptide);
    
    // Create PeptideHit with the same sequence
    PeptideHit hit = AScoreTestHelper::generatePeptideHit("STYSTSYSTSTYS(Phospho)TSYT(Phospho)");
    
    // Run AScore
    hit = ascore.compute(hit, spectrum);
    
    // Verify results
    TEST_EQUAL(hit.getSequence().toString(), "STYSTSYSTSTYS(Phospho)TSYT(Phospho)");
    TEST_EQUAL(hit.metaValueExists("ProForma"), true);
    
    // Verify ProForma string contains the phosphorylated sites
    String proforma = hit.getMetaValue("ProForma");
    TEST_EQUAL(proforma, "STYSTSYSTSTYS[Phospho|score=0.9780]TSYT[Phospho|score=0.9780]");

  }
  
  // Test case 8: Peptide with no phosphorylation sites
  {
    // Create a peptide with no S/T/Y residues
    AASequence peptide = AASequence::fromString("PEPNIDEA");
    
    // Generate theoretical spectrum
    PeakSpectrum spectrum = AScoreTestHelper::generateTheoreticalSpectrum(peptide);
    
    // Create PeptideHit with the same sequence
    PeptideHit hit = AScoreTestHelper::generatePeptideHit("PEPNIDEA");
    
    // Run AScore
    hit = ascore.compute(hit, spectrum);
    
    // Verify results - expect score of -1 for no phosphorylation
    TEST_EQUAL(hit.getScore(), -1.0);
    
    // ProForma should not exist since there are no phosphorylation sites
    TEST_EQUAL(hit.metaValueExists("ProForma"), false);
  }
  
  // Test case 9: Peptide with a large number of permutations
  {
    // Create a peptide with many S/T/Y residues and multiple phosphorylations
    // This will generate a large number of permutations
    AASequence peptide = AASequence::fromString("STYSTSY(Phospho)STSTYS(Phospho)TSYT(Phospho)");
    
    // Set max_num_perm to a low value to test early termination
    params.setValue("max_num_perm", 10);
    ascore.setParameters(params);
    
    // Generate theoretical spectrum
    PeakSpectrum spectrum = AScoreTestHelper::generateTheoreticalSpectrum(peptide);
    
    // Create PeptideHit with the same sequence
    PeptideHit hit = AScoreTestHelper::generatePeptideHit("STYSTSY(Phospho)STSTYS(Phospho)TSYT(Phospho)");
    
    // Run AScore
    hit = ascore.compute(hit, spectrum);
    
    // Reset max_num_perm for other tests
    params.setValue("max_num_perm", 16384);
    ascore.setParameters(params);
  }
  
  // Test case 10: PhosphoDecoy sites with PhosphoDecoy disabled
  {
    // Disable PhosphoDecoy
    params.setValue("add_decoys", "false");
    ascore.setParameters(params);
    
    // Create peptide with PhosphoDecoy site
    AASequence peptide = AASequence::fromString("PEPTIDEA(PhosphoDecoy)");
    
    // Generate theoretical spectrum
    PeakSpectrum spectrum = AScoreTestHelper::generateTheoreticalSpectrum(peptide);
    
    // Create PeptideHit with the same sequence
    PeptideHit hit = AScoreTestHelper::generatePeptideHit("PEPTIDEA(PhosphoDecoy)");
    
    // Run AScore - should return the original hit without processing
    hit = ascore.compute(hit, spectrum);
    
    // Verify that the hit was not processed (score should be -1)
    TEST_EQUAL(hit.getScore(), -1.0);
  }
  
  // Reset parameters for other tests
  params.setValue("add_decoys", "false");
  ascore.setParameters(params);
}
END_SECTION

delete ptr_test;
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
