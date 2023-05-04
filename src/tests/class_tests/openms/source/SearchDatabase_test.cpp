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
// $Maintainer: $
// $Authors: Max Alcer, Heike Einsfeld $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/ID/SearchDatabase.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <vector>
#include <OpenMS/MATH/MISC/MathFunctions.h>

using namespace OpenMS;
using namespace std;

class SearchDatabase_test : public SearchDatabase
{
  public:

  using SearchDatabase::SearchDatabase;

  void printAllFragments()
  {
  
  int counter = 0;

  for (auto j : all_fragments_)
  { 
    if(counter == bucketsize_)
    {
      cout << "\n";
      counter = 0;
    }

    cout << j.fragment_mz_ << " " << all_peptides_[j.peptide_index_].peptide_mz_ << "\n";

    counter++;
  }

  }
  
  std::vector<Fragment_> getAllFragments(){return all_fragments_;}

  size_t getBucketsize(){return bucketsize_;}

  string getDigestorEnzyme(){return digestor_enzyme_;}

  size_t getMissedCleavages(){return missed_cleavages_;}
  	
  size_t getpeptide_minlen(){return peptide_min_length_;}
  size_t getpeptide_maxlen(){return peptide_max_length_;}
  double getpeptide_minMass(){return peptide_min_mass_;}
  double getpeptide_maxMass(){return peptide_max_mass_;}
  double getfrag_minMZ(){return fragment_min_mz_;}
  double getfrag_maxMZ(){return fragment_max_mz_;}
  size_t getpeptideSize(){return all_peptides_.size();}
  vector<SearchDatabase::Peptide_> getpeptide(){return all_peptides_;}

  double getFragMZ(size_t i)
  {
    return all_fragments_[i].fragment_mz_;
  }

  double getPrecMZ(size_t i)
  {
    return all_peptides_[all_fragments_[i].peptide_index_].peptide_mz_;
  }

};

START_TEST(SearchDatabase, "$Id:$")


START_SECTION(SearchDatabase(const std::vector<FASTAFile::FASTAEntry>& entries))

  std::vector<FASTAFile::FASTAEntry> entries{
  {"test1", "test1", "LRLRACGLNFADLMARQGLY"},
  {"test2", "test2", "AAASPPLLRCLVLTGFGGYD"},
  {"test3", "test3", "KVKLQSRPAAPPAPGPGQLT"}};

  SearchDatabase_test sdb(entries);

  vector<AASequence> all_peptides;

  ProteaseDigestion digestor;
  
  digestor.setEnzyme(sdb.getDigestorEnzyme());
  digestor.setMissedCleavages(sdb.getMissedCleavages());

  for (auto i : entries)
  {
    vector<AASequence> peptides;

    digestor.digest(AASequence::fromString(i.sequence), peptides, sdb.getpeptide_minlen(), sdb.getpeptide_maxlen());
    for (const auto& pep : peptides)
    { 
      if (pep.toString().find('X') != string::npos) continue;
      double seq_mz = pep.getMonoWeight();
      if (!Math::contains(seq_mz, sdb.getpeptide_minMass(), sdb.getpeptide_maxMass())) continue;
      all_peptides.emplace_back(pep);
        
    }
  }

  TEST_EQUAL(all_peptides.size(), sdb.getpeptideSize())
  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum b_y_ions;
  int count_all_frags=0;  
    
  for (size_t i = 0; i < all_peptides.size(); i++)
  { 
    tsg.getSpectrum(b_y_ions, all_peptides[i], 1, 1);      
    for (const auto& frag : b_y_ions)
    { 
      if (!Math::contains(frag.getMZ(), sdb.getfrag_minMZ(), sdb.getfrag_maxMZ())) continue;
      count_all_frags++;        
    }
    b_y_ions.clear(true);
  }
  TEST_EQUAL(sdb.getAllFragments().size(), count_all_frags)
  
END_SECTION

START_SECTION(void search(MSSpectrum& spectrum, std::vector<Candidate>& candidates))  

  cout << "\n";

  std::vector<FASTAFile::FASTAEntry> entries{
  {"test1", "test1", "LRLRACGLNFADLMARQGLY"},
  {"test2", "test2", "AAASPPLLRCLVLTGFGGYD"},
  {"test3", "test3", "KVKLQSRPAAPPAPGPGQLT"}};

  SearchDatabase_test sdb(entries);

  MSSpectrum spec;

  Precursor prec{};

  std::vector<SearchDatabase::Candidate> candidates;

  START_SECTION(Searching 3 Fragments it should find)

  prec.setCharge(1);

  prec.setMZ(1280.6);

  spec.setPrecursors({prec});

  spec.push_back({605.308, 100});
  spec.push_back({676.345, 100});
  spec.push_back({823.413, 100});  

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 1)

  END_SECTION

  START_SECTION(Searching Fragment it should not find because of Fragment Mass)

  spec.clear(true);

  spec.setPrecursors({prec});

  spec.push_back({1040, 100});

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

  START_SECTION(Searching Fragment it should not find because of Precursor Mass)

  spec.clear(true);

  prec.setMZ(1500);

  spec.setPrecursors({prec});

  spec.push_back({572.304, 100});

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

  START_SECTION(Searching Fragment it should not find because its smaller then all Fragments in Database)

  spec.clear(true);

  prec.setMZ(1500);

  spec.setPrecursors({prec});

  spec.push_back({100, 100});

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

  START_SECTION(Searching Fragment it should not find because its bigger then all Fragments in Database)

  spec.clear(true);

  prec.setMZ(1500);

  spec.setPrecursors({prec});

  spec.push_back({2000, 100});

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

  START_SECTION(Testing filtering of Spectrum by best Peaks)

  spec.clear(true);

  prec.setMZ(1280.6);

  spec.setPrecursors({prec});

  spec.push_back({2000, 80});

  spec.push_back({1500, 99});

  spec.push_back({500, 85});

  spec.push_back({1000, 90});

  spec.push_back({100, 100});

  spec.push_back({937.456, 5});

  sdb.search(spec, candidates);

  TEST_EQUAL(candidates.size(), 0)

  END_SECTION

  // sdb.printAllFragments();

END_SECTION

START_SECTION(void search(MSExperiment& experiment, std::vector<std::pair<std::vector<Candidate>, size_t>>& candidates))

  std::vector<FASTAFile::FASTAEntry> entries{
  {"test1", "test1", "LRLRACGLNFADLMARQGLY"},
  {"test2", "test2", "AAASPPLLRCLVLTGFGGYD"},
  {"test3", "test3", "KVKLQSRPAAPPAPGPGQLT"}};

  SearchDatabase_test sdb(entries);

  MSExperiment exp;

  MSSpectrum spec;

  Precursor prec{};

  prec.setCharge(1);

  prec.setMZ(1281.6);

  spec.setPrecursors({prec});

  spec.push_back({605.318, 100});

  exp.addSpectrum(spec);

  spec.clear(true);

  prec.setMZ(894.529);

  spec.setPrecursors({prec});

  spec.push_back({175.119, 100});

  exp.addSpectrum(spec);

  spec.clear(true);

  prec.setMZ(1655.89);

  spec.setPrecursors({prec});

  spec.push_back({1544.83, 100});

  exp.addSpectrum(spec);  

  std::vector<std::pair<std::vector<SearchDatabase::Candidate>, size_t>> candidates;

  sdb.search(exp, candidates);

  TEST_EQUAL(candidates.size(), 3)

END_SECTION

END_TEST