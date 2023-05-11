// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Chris Bielow $
// $Authors: Max Alcer, Heike Einsfeld $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <vector>

using namespace OpenMS;
using namespace Constants::UserParam;
using namespace std;

class FragmentIndex_test : public FragmentIndex
{
  public:

  using FragmentIndex::FragmentIndex;
  
  std::vector<Fragment_> getAllFragments(){return all_fragments_;}

  bool isSortedBucketFragsMZ()
  {
    return is_sorted(bucket_frags_mz_.begin(), bucket_frags_mz_.end());
  }

  bool isSortedAllFragments()
  {
    bool test_sorted = true;

    for (size_t i = 0; i < all_fragments_.size(); i += bucketsize_)
    {
      auto bucket_begin = all_fragments_.begin()+i;
      auto condition = distance(all_fragments_.begin(), bucket_begin+bucketsize_) >= int(all_fragments_.size());

      test_sorted &= is_sorted(bucket_begin, condition ? all_fragments_.end() : bucket_begin+bucketsize_, 
      [&](const Fragment_& l, const Fragment_& r) -> bool 
      {return (all_peptides_[l.peptide_index_].peptide_mz_ < all_peptides_[r.peptide_index_].peptide_mz_);});
    }
    return test_sorted;
  }

};

START_TEST(FragmentIndex, "$Id:$")


START_SECTION(build(const std::vector<FASTAFile::FASTAEntry>& entries))

  std::vector<FASTAFile::FASTAEntry> entries{
  {"test1", "test1", "LRLRACGLNFADLMARQGLY"},
  {"test2", "test2", "AAASPPLLRCLVLTGFGGYD"},
  {"test3", "test3", "KVKLQSRPAAPPAPGPGQLT"}};

  FragmentIndex_test sdb;
  
  START_SECTION(test number of fragments)
    sdb.build(entries);
    TEST_EQUAL(187, sdb.getAllFragments().size())
  END_SECTION

  START_SECTION(test sortation)
    TEST_TRUE(sdb.isSortedBucketFragsMZ())

    TEST_TRUE(sdb.isSortedAllFragments())
  END_SECTION

  START_SECTION(test without modifications)

  auto params = sdb.getParameters();

  params.setValue("modifications_fixed", std::vector<std::string>{});
  params.setValue("modifications_variable", std::vector<std::string>{});

  sdb.setParameters(params);

  sdb.build(entries);
  TEST_NOT_EQUAL(187, sdb.getAllFragments().size());

  END_SECTION

  START_SECTION(test peptide mz bounds)

  auto params = sdb.getParameters();

  params.setValue("peptide_min_mass", 2000);
  params.setValue("peptide_max_mass", 2000);

  sdb.setParameters(params);

  sdb.build(entries);
  TEST_EQUAL(0, sdb.getAllFragments().size());

  END_SECTION

  START_SECTION(test fragment mz bounds)

  auto params = sdb.getParameters();

  params.setValue("fragment_min_mz", 500);
  params.setValue("fragment_max_mz", 500);

  sdb.setParameters(params);

  sdb.build(entries);
  TEST_EQUAL(0, sdb.getAllFragments().size());

  END_SECTION
  
END_SECTION

START_SECTION(void search(MSSpectrum& spectrum, std::vector<Candidate>& candidates))  

  cout << "\n";

  std::vector<FASTAFile::FASTAEntry> entries{
  {"test1", "test1", "LRLRACGLNFADLMARQGLY"},
  {"test2", "test2", "AAASPPLLRCLVLTGFGGYD"},
  {"test3", "test3", "KVKLQSRPAAPPAPGPGQLT"}};

  FragmentIndex_test sdb;
  std::vector<FragmentIndex::Candidate> candidates;

  MSSpectrum spec;

  Precursor prec{};

  START_SECTION(searching without building)

  prec.setCharge(1);

  prec.setMZ(1281.6);

  spec.setPrecursors({prec});
  spec.setMSLevel(2);

  spec.push_back({605.308, 100});
  spec.push_back({676.345, 100});
  spec.push_back({823.413, 100});  

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

  sdb.build(entries);

  START_SECTION(Searching 3 Fragments it should find (with Da and ppm))

  prec.setCharge(1);

  prec.setMZ(1281.6);

  spec.setPrecursors({prec});
  spec.setMSLevel(2);

  spec.push_back({605.308, 100});
  spec.push_back({676.345, 100});
  spec.push_back({823.413, 100});  

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 1)

  auto params = sdb.getParameters();

  params.setValue("fragment_mz_tolerance_unit", UNIT_PPM);
  params.setValue("fragment_mz_tolerance", 5.f);
  params.setValue("precursor_mz_tolerance_unit", UNIT_PPM);
  params.setValue("precursor_mz_tolerance", 50.f);

  sdb.setParameters(params);

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 1)

  END_SECTION

  START_SECTION(Searching Fragment it should not find because of Fragment Mass)

  auto params = sdb.getParameters();

  params.setValue("fragment_mz_tolerance_unit", UNIT_DA);
  params.setValue("fragment_mz_tolerance", 0.05);
  params.setValue("precursor_mz_tolerance_unit", UNIT_DA);
  params.setValue("precursor_mz_tolerance", 2.0);

  sdb.setParameters(params);

  spec.clear(false);

  spec.push_back({1040, 100});

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

  START_SECTION(Searching Fragment it should not find because of Precursor Mass)

  spec.clear(true);

  prec.setMZ(1500);

  spec.setPrecursors({prec});
  spec.setMSLevel(2);

  spec.push_back({572.304, 100});

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

  START_SECTION(Searching Fragment it should not find because its smaller then all Fragments in Database)

  spec.clear(false);

  spec.push_back({100, 100});

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

  START_SECTION(Searching Fragment it should not find because its bigger then all Fragments in Database)

  spec.clear(false);

  spec.push_back({2000, 100});

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

  START_SECTION(Testing filtering of Spectrum by best Peaks)

  spec.clear(true);

  prec.setMZ(1281.6);

  spec.setPrecursors({prec});
  spec.setMSLevel(2);

  spec.push_back({2000, 80});

  spec.push_back({1500, 99});

  spec.push_back({500, 85});

  spec.push_back({1000, 90});

  spec.push_back({100, 100});

  spec.push_back({937.456, 5});

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

END_SECTION

START_SECTION(void search(MSExperiment& experiment, std::vector<CandidatesWithIndex>& candidates))

  std::vector<FASTAFile::FASTAEntry> entries{
  {"test1", "test1", "LRLRACGLNFADLMARQGLY"},
  {"test2", "test2", "AAASPPLLRCLVLTGFGGYD"},
  {"test3", "test3", "KVKLQSRPAAPPAPGPGQLT"}};

  FragmentIndex_test sdb;

  sdb.build(entries);

  MSExperiment exp;

  MSSpectrum spec;

  Precursor prec{};

  prec.setCharge(1);

  prec.setMZ(1281.6);

  spec.setPrecursors({prec});
  spec.setMSLevel(2);

  spec.push_back({605.318, 100});

  exp.addSpectrum(spec);

  spec.clear(true);

  prec.setMZ(894.529);

  spec.setPrecursors({prec});
  spec.setMSLevel(2);

  spec.push_back({175.119, 100});

  exp.addSpectrum(spec);

  spec.clear(true);

  prec.setMZ(1655.89);

  spec.setPrecursors({prec});
  spec.setMSLevel(2);

  spec.push_back({1544.83, 100});

  exp.addSpectrum(spec);  

  std::vector<FragmentIndex::CandidatesWithIndex> candidates;

  sdb.search(exp, candidates);

  TEST_EQUAL(candidates.size(), 3)

END_SECTION

END_TEST