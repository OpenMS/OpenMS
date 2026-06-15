// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/HydrophobicityProfile.h>

///////////////////////////

#include <OpenMS/CONCEPT/Constants.h>
#include <cmath>
#include <map>

using namespace OpenMS;
using namespace std;

START_TEST(HydrophobicityProfile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HydrophobicityProfile* ptr = nullptr;
HydrophobicityProfile* null_ptr = nullptr;
START_SECTION(HydrophobicityProfile())
{
  ptr = new HydrophobicityProfile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~HydrophobicityProfile())
{
  delete ptr;
}
END_SECTION

START_SECTION(double computeGRAVY(const AASequence& seq))
{
  AASequence seq("ACDE");
  HydrophobicityProfile profile;
  double gravy = profile.computeGRAVY(seq);
  TEST_REAL_SIMILAR(gravy, -0.675); 
  AASequence seq_2;
  TEST_EXCEPTION(Exception::InvalidValue,profile.computeGRAVY(seq_2)); 
  AASequence seq_3("XXX");
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue,profile.computeGRAVY(seq_3),
  "the value 'X' was used but is not valid; No hydrophobicity value known for this residue");
}
END_SECTION

START_SECTION(std::vector<double> computeProfile(const AASequence& seq, const HydrophobicityScaleMethod scale) const)
{
  HydrophobicityProfile profile;
  AASequence seq("ACDE");
  AASequence seq_2;
  AASequence seq_3("XXX");
  std::vector<double> vec = profile.computeProfile(seq,HydrophobicityScaleMethod::KYTE_DOOLITTLE);
  TEST_REAL_SIMILAR(vec[0],1.8);
  TEST_REAL_SIMILAR(vec[1],2.5);
  TEST_REAL_SIMILAR(vec[2],-3.5);
  TEST_REAL_SIMILAR(vec[3],-3.5); 
  TEST_REAL_SIMILAR(profile.computeProfile(seq,HydrophobicityScaleMethod::EISENBERG)[0],0.62);
  TEST_EXCEPTION(Exception::InvalidValue,profile.computeProfile(seq_3,HydrophobicityScaleMethod::EISENBERG))
}
END_SECTION

START_SECTION(std::vector<double> computeWindowedProfile(const AASequence& seq, Size window_size, const HydrophobicityScaleMethod scale) const)
{
  HydrophobicityProfile profile;
  AASequence seq("ACDEF");
  AASequence seq_2;
  std::vector<double> vec = profile.computeWindowedProfile(seq,3);
  TEST_REAL_SIMILAR(vec[0],0.266666666666667);
  TEST_REAL_SIMILAR(vec[1],-1.5);
  TEST_REAL_SIMILAR(vec[2],-1.4);
  TEST_REAL_SIMILAR(profile.computeWindowedProfile(seq,6)[0],0.02);
  TEST_EXCEPTION(Exception::InvalidSize,profile.computeWindowedProfile(seq,0));
  TEST_EXCEPTION(Exception::InvalidValue,profile.computeWindowedProfile(seq_2,3));
}
END_SECTION

START_SECTION(std::vector<double> computeHydrophobicMoment(const AASequence& seq, Size window_size, double angle) const)
{
  HydrophobicityProfile profile;
  AASequence seq_1("ACDEF");
  AASequence seq_2;
  std::vector<double> vec = profile.computeHydrophobicMoment(seq_1,3,100);
  TEST_REAL_SIMILAR(vec[0],0.511576803);
  TEST_REAL_SIMILAR(vec[1],0.435170599);
  TEST_REAL_SIMILAR(vec[2],0.734926405);
  TEST_EXCEPTION(Exception::InvalidSize,profile.computeHydrophobicMoment(seq_1,0));
  TEST_EXCEPTION(Exception::InvalidValue,profile.computeHydrophobicMoment(seq_2,3));
}
END_SECTION

START_SECTION([EXTRA] computeHydrophobicMoment matches the published Eisenberg reference)
{
  // Cross-check computeHydrophobicMoment against an independent computation of the
  // Eisenberg (1984) hydrophobic moment from the *published* Eisenberg consensus
  // hydrophobicity scale (cf. R package R.Peptides, hmoment.R). The section above
  // pins OpenMS's own output; this validates that output against the literature
  // scale + formula:  moment_w = sqrt( (sum h_i sin(i*d))^2 + (sum h_i cos(i*d))^2 ) / w
  HydrophobicityProfile profile;

  // Published Eisenberg consensus hydrophobicity values (Eisenberg et al. 1984).
  const std::map<char, double> eisenberg = {
    {'A',  0.62}, {'R', -2.53}, {'N', -0.78}, {'D', -0.90}, {'C',  0.29},
    {'Q', -0.85}, {'E', -0.74}, {'G',  0.48}, {'H', -0.40}, {'I',  1.38},
    {'L',  1.06}, {'K', -1.50}, {'M',  0.64}, {'F',  1.19}, {'P',  0.12},
    {'S', -0.18}, {'T', -0.05}, {'W',  0.81}, {'Y',  0.26}, {'V',  1.08}};

  auto reference_moment = [&](const std::string& seq, Size window, double angle_deg)
  {
    const double a = angle_deg * Constants::PI / 180.0;
    std::vector<double> out;
    for (Size i = 0; i + window <= seq.size(); ++i)
    {
      double sum_sin = 0.0, sum_cos = 0.0;
      for (Size pos = 0; pos < window; ++pos)
      {
        const double h = eisenberg.at(seq[i + pos]);
        sum_sin += h * std::sin(a * pos);
        sum_cos += h * std::cos(a * pos);
      }
      out.push_back(std::sqrt(sum_sin * sum_sin + sum_cos * sum_cos) / window);
    }
    return out;
  };

  struct Case { std::string pep; Size window; double angle; };
  const std::vector<Case> cases = {
    {"ACDEF",     3, 100.0},   // also the peptide pinned literally in the section above
    {"FLIGKWVRP", 3, 100.0},
    {"FLIGKWVRP", 5, 160.0}};

  for (const Case& c : cases)
  {
    std::vector<double> got = profile.computeHydrophobicMoment(AASequence::fromString(c.pep), c.window, c.angle);
    std::vector<double> ref = reference_moment(c.pep, c.window, c.angle);
    TEST_EQUAL(got.size(), ref.size())
    for (Size i = 0; i < ref.size() && i < got.size(); ++i)
    {
      TEST_REAL_SIMILAR(got[i], ref[i])
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST