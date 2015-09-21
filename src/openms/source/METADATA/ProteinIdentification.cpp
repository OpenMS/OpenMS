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
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <sstream>
#include <algorithm>

using namespace std;

namespace OpenMS
{

  const std::string ProteinIdentification::NamesOfPeakMassType[] = {"Monoisotopic", "Average"};
  const std::string ProteinIdentification::NamesOfDigestionEnzyme[] = {"Trypsin", "Pepsin A", "Protease K", "Chymotrypsin", "No enzyme", "Unknown"};

  ProteinIdentification::ProteinGroup::ProteinGroup() :
    probability(0.0), accessions()
  {
  }

  bool ProteinIdentification::ProteinGroup::operator==(const ProteinGroup& rhs) const
  {
    return (probability == rhs.probability) && (accessions == rhs.accessions);
  }

  bool ProteinIdentification::ProteinGroup::operator<(const ProteinGroup& rhs) const
  {
    // comparison of probabilities is intentionally "the wrong way around":
    if (probability > rhs.probability) return true;
    if (probability < rhs.probability) return false;
    if (accessions.size() < rhs.accessions.size()) return true;
    if (accessions.size() > rhs.accessions.size()) return false;
    return accessions < rhs.accessions;
  }


  ProteinIdentification::SearchParameters::SearchParameters() :
    db(),
    db_version(),
    taxonomy(),
    charges(),
    mass_type(MONOISOTOPIC),
    fixed_modifications(),
    variable_modifications(),
    enzyme(UNKNOWN_ENZYME),
    missed_cleavages(0),
    peak_mass_tolerance(0.0),
    precursor_tolerance(0.0),
    digestion_enzyme("unknown_enzyme","")
  {
  }

  bool ProteinIdentification::SearchParameters::operator==(const SearchParameters& rhs) const
  {
    return db == rhs.db &&
           db_version == rhs.db_version &&
           taxonomy == rhs.taxonomy &&
           charges == rhs.charges &&
           mass_type == rhs.mass_type &&
           fixed_modifications == rhs.fixed_modifications &&
           variable_modifications == rhs.variable_modifications &&
           enzyme == rhs.enzyme &&
           missed_cleavages == rhs.missed_cleavages &&
           peak_mass_tolerance == rhs.peak_mass_tolerance &&
           precursor_tolerance == rhs.precursor_tolerance &&
           digestion_enzyme == rhs.digestion_enzyme;
  }

  bool ProteinIdentification::SearchParameters::operator!=(const SearchParameters& rhs) const
  {
    return !(*this == rhs);
  }

  ProteinIdentification::ProteinIdentification() :
    MetaInfoInterface(),
    id_(),
    search_engine_(),
    search_engine_version_(),
    search_parameters_(),
    date_(),
    protein_score_type_(),
    higher_score_better_(true),
    protein_hits_(),
    protein_groups_(),
    indistinguishable_proteins_(),
    protein_significance_threshold_(0.0)
  {
  }

  ProteinIdentification::ProteinIdentification(const ProteinIdentification& source) :
    MetaInfoInterface(source),
    id_(source.id_),
    search_engine_(source.search_engine_),
    search_engine_version_(source.search_engine_version_),
    search_parameters_(source.search_parameters_),
    date_(source.date_),
    protein_score_type_(source.protein_score_type_),
    higher_score_better_(source.higher_score_better_),
    protein_hits_(source.protein_hits_),
    protein_groups_(source.protein_groups_),
    indistinguishable_proteins_(source.indistinguishable_proteins_),
    protein_significance_threshold_(source.protein_significance_threshold_)
  {
  }

  ProteinIdentification::~ProteinIdentification()
  {
  }

  void ProteinIdentification::setDateTime(const DateTime& date)
  {
    date_ = date;
  }

  const DateTime& ProteinIdentification::getDateTime() const
  {
    return date_;
  }

  const vector<ProteinHit>& ProteinIdentification::getHits() const
  {
    return protein_hits_;
  }

  vector<ProteinHit>& ProteinIdentification::getHits()
  {
    return protein_hits_;
  }

  void ProteinIdentification::setHits(const vector<ProteinHit>& protein_hits)
  {
    protein_hits_ = protein_hits;
  }

  vector<ProteinHit>::iterator ProteinIdentification::findHit(
    const String& accession)
  {
    vector<ProteinHit>::iterator pos = protein_hits_.begin();
    for (; pos != protein_hits_.end(); ++pos)
    {
      if (pos->getAccession() == accession)
        break;
    }
    return pos;
  }

  const vector<ProteinIdentification::ProteinGroup>& ProteinIdentification::getProteinGroups() const
  {
    return protein_groups_;
  }

  vector<ProteinIdentification::ProteinGroup>& ProteinIdentification::getProteinGroups()
  {
    return protein_groups_;
  }

  void ProteinIdentification::insertProteinGroup(const ProteinIdentification::ProteinGroup& group)
  {
    protein_groups_.push_back(group);
  }

  const vector<ProteinIdentification::ProteinGroup>&
  ProteinIdentification::getIndistinguishableProteins() const
  {
    return indistinguishable_proteins_;
  }

  vector<ProteinIdentification::ProteinGroup>&
  ProteinIdentification::getIndistinguishableProteins()
  {
    return indistinguishable_proteins_;
  }

  void ProteinIdentification::insertIndistinguishableProteins(
    const ProteinIdentification::ProteinGroup& group)
  {
    indistinguishable_proteins_.push_back(group);
  }

  // retrieval of the peptide significance threshold value
  double ProteinIdentification::getSignificanceThreshold() const
  {
    return protein_significance_threshold_;
  }

  // setting of the peptide significance threshold value
  void ProteinIdentification::setSignificanceThreshold(double value)
  {
    protein_significance_threshold_ = value;
  }

  void ProteinIdentification::setScoreType(const String& type)
  {
    protein_score_type_ = type;
  }

  const String& ProteinIdentification::getScoreType() const
  {
    return protein_score_type_;
  }

  void ProteinIdentification::insertHit(const ProteinHit& protein_hit)
  {
    protein_hits_.push_back(protein_hit);
  }

  ProteinIdentification& ProteinIdentification::operator=(const ProteinIdentification& source)
  {
    if (this == &source)
    {
      return *this;
    }
    MetaInfoInterface::operator=(source);
    id_ = source.id_;
    search_engine_ = source.search_engine_;
    search_engine_version_ = source.search_engine_version_;
    search_parameters_ = source.search_parameters_;
    date_ = source.date_;
    protein_hits_ = source.protein_hits_;
    protein_groups_ = source.protein_groups_;
    indistinguishable_proteins_ = source.indistinguishable_proteins_;
    protein_score_type_ = source.protein_score_type_;
    protein_significance_threshold_ = source.protein_significance_threshold_;
    higher_score_better_ = source.higher_score_better_;
    return *this;
  }

  // Equality operator
  bool ProteinIdentification::operator==(const ProteinIdentification& rhs) const
  {
    return MetaInfoInterface::operator==(rhs) &&
           id_ == rhs.id_ &&
           search_engine_ == rhs.search_engine_ &&
           search_engine_version_ == rhs.search_engine_version_ &&
           search_parameters_ == rhs.search_parameters_ &&
           date_ == rhs.date_ &&
           protein_hits_ == rhs.protein_hits_ &&
           protein_groups_ == rhs.protein_groups_ &&
           indistinguishable_proteins_ == rhs.indistinguishable_proteins_ &&
           protein_score_type_ == rhs.protein_score_type_ &&
           protein_significance_threshold_ == rhs.protein_significance_threshold_ &&
           higher_score_better_ == rhs.higher_score_better_;

  }

  // Inequality operator
  bool ProteinIdentification::operator!=(const ProteinIdentification& rhs) const
  {
    return !operator==(rhs);
  }

  void ProteinIdentification::sort()
  {
    if (higher_score_better_)
    {
      std::sort(protein_hits_.begin(), protein_hits_.end(), ProteinHit::ScoreMore());
    }
    else
    {
      std::sort(protein_hits_.begin(), protein_hits_.end(), ProteinHit::ScoreLess());
    }
  }

  void ProteinIdentification::assignRanks()
  {
    if (protein_hits_.empty())
      return;

    UInt rank = 1;
    sort();
    vector<ProteinHit>::iterator lit = protein_hits_.begin();
    float tmpscore = lit->getScore();
    while (lit != protein_hits_.end())
    {
      lit->setRank(rank);
      ++lit;
      if (lit != protein_hits_.end() && lit->getScore() != tmpscore)
      {
        ++rank;
        tmpscore = lit->getScore();
      }
    }
  }

  Size ProteinIdentification::computeCoverage(const std::vector<PeptideIdentification>& pep_ids)
  {
    // TODO: we currently ignore overlapping peptides, i.e. the coverage could be > 100%

    Size no_seq_count(0);
    // index the proteins by accession
    // Accession -> set of pep sequences
    // (use set to discard multi-pep matches)
    Map<String, std::set<String> > protein_index;
    for (Size i = 0; i < protein_hits_.size(); ++i)
    {
      std::set<String> empty;
      protein_index[protein_hits_[i].getAccession()] = empty;
      if (protein_hits_[i].getSequence().length() == 0)
        ++no_seq_count;
    }

    if (no_seq_count > 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, String(no_seq_count) + " of " + protein_hits_.size() + " ProteinHits do not contain a protein sequence. Cannot compute coverage! Use PeptideIndexer to annotate proteins with sequence information.");
    }

    // go through peptides and add length to proteinHit
    Size protein_not_found_counter(0);
    for (vector<PeptideIdentification>::const_iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
    {
      // peptide hits
      vector<PeptideHit> peptide_hits = it1->getHits();
      for (vector<PeptideHit>::iterator it2 = peptide_hits.begin(); it2 != peptide_hits.end(); ++it2)
      {
        set<String> protein_accessions = it2->extractProteinAccessions();
        // matched proteins for hit
        for (set<String>::const_iterator it3 = protein_accessions.begin(); it3 != protein_accessions.end(); ++it3)
        {
          if (protein_index.has(*it3))
          {
            protein_index[*it3].insert(it2->getSequence().toUnmodifiedString());
          }
          else
          {
            ++protein_not_found_counter;
          }
        }
      }
    }

    if (protein_not_found_counter > 0)
      LOG_WARN << "ProteinIdentification::computeCoverage() was given PeptideIdentifications where " << protein_not_found_counter << " did not match a known Protein!" << std::endl;

    // store coverage
    for (Size i = 0; i < protein_hits_.size(); ++i)
    {
      // add up peptide sizes
      Size covered_length(0);
      for (std::set<String>::const_iterator it = protein_index[protein_hits_[i].getAccession()].begin();
           it != protein_index[protein_hits_[i].getAccession()].end();
           ++it)
      {
        covered_length += it->size();
      }
      // set coverage
      protein_hits_[i].setCoverage(double(covered_length) / (double)protein_hits_[i].getSequence().length() * 100.0);
    }

    return protein_not_found_counter;
  }

  bool ProteinIdentification::isHigherScoreBetter() const
  {
    return higher_score_better_;
  }

  void ProteinIdentification::setHigherScoreBetter(bool value)
  {
    higher_score_better_ = value;
  }

  const String& ProteinIdentification::getIdentifier() const
  {
    return id_;
  }

  void ProteinIdentification::setIdentifier(const String& id)
  {
    id_ = id;
  }

  void ProteinIdentification::setSearchEngine(const String& search_engine)
  {
    search_engine_ = search_engine;
  }

  const String& ProteinIdentification::getSearchEngine() const
  {
    return search_engine_;
  }

  void ProteinIdentification::setSearchEngineVersion(const String& search_engine_version)
  {
    search_engine_version_ = search_engine_version;
  }

  const String& ProteinIdentification::getSearchEngineVersion() const
  {
    return search_engine_version_;
  }

  void ProteinIdentification::setSearchParameters(const SearchParameters& search_parameters)
  {
    search_parameters_ = search_parameters;
  }

  const ProteinIdentification::SearchParameters& ProteinIdentification::getSearchParameters() const
  {
    return search_parameters_;
  }

} // namespace OpenMS
