// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ProteinHit.h>
#include <ostream>


using namespace std;

namespace OpenMS
{
  const double ProteinHit::COVERAGE_UNKNOWN = -1;

  // default constructor
  ProteinHit::ProteinHit() :
    MetaInfoInterface(),
    score_(0),
    rank_(0),
    accession_(""),
    sequence_(""),
    coverage_(COVERAGE_UNKNOWN)
  {
  }

  // values constructor
  ProteinHit::ProteinHit(double score, UInt rank, String accession, String sequence) :
    MetaInfoInterface(),
    score_(score),
    rank_(rank),
    accession_(accession.trim()),
    sequence_(sequence.trim()),
    coverage_(COVERAGE_UNKNOWN)
  {
  }

  // assignment operator for MetaInfoInterface
  ProteinHit& ProteinHit::operator=(const MetaInfoInterface& source)
  {
    MetaInfoInterface::operator=(source);
    return *this;
  }

  // equality operator
  bool ProteinHit::operator==(const ProteinHit& rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && score_ == rhs.score_
           && rank_ == rhs.rank_
           && accession_ == rhs.accession_
           && sequence_ == rhs.sequence_
           && coverage_ == rhs.coverage_
           && modifications_ == rhs.modifications_;
  }

  // inequality operator
  bool ProteinHit::operator!=(const ProteinHit& rhs) const
  {
    return !operator==(rhs);
  }

  // returns the score of the protein hit
  double ProteinHit::getScore() const
  {
    return score_;
  }

  // returns the rank of the protein hit
  UInt ProteinHit::getRank() const
  {
    return rank_;
  }

  // returns the protein sequence
  const String& ProteinHit::getSequence() const
  {
    return sequence_;
  }

  // returns the accession of the protein
  const String& ProteinHit::getAccession() const
  {
    return accession_;
  }

  // returns the description of the protein
  String ProteinHit::getDescription() const
  {
    return getMetaValue("Description").toString();
  }

  // returns the coverage (in percent) of the protein hit based upon matched peptides
  double ProteinHit::getCoverage() const
  {
    return coverage_;
  }

  // sets the score of the protein hit
  void ProteinHit::setScore(const double score)
  {
    score_ = score;
  }

  // sets the rank
  void ProteinHit::setRank(UInt newrank)
  {
    rank_ = newrank;
  }

  // sets the protein sequence
  void ProteinHit::setSequence(const String& sequence)
  {
    sequence_ = sequence;
    sequence_.trim();
  }

  // sets the protein sequence
  void ProteinHit::setSequence(String&& sequence)
  {
    sequence_ = std::move(sequence);
    sequence_.trim();
  }


  // sets the description of the protein
  void ProteinHit::setDescription(const String& description)
  {
    setMetaValue("Description", description);
  }

  // sets the accession of the protein
  void ProteinHit::setAccession(const String& accession)
  {
    accession_ = accession;
    accession_.trim();
  }

  // sets the coverage (in percent) of the protein hit based upon matched peptides
  void ProteinHit::setCoverage(const double coverage)
  {
    coverage_ = coverage;
  }

  const set<pair<Size, ResidueModification>>& ProteinHit::getModifications() const
  {
    return modifications_;
  }

  void ProteinHit::setModifications(std::set<std::pair<Size, ResidueModification>>& mods)
  {
    modifications_ = mods;
  }

  std::ostream& operator<< (std::ostream& stream, const ProteinHit& hit)
  {
    return stream << "protein hit with accession '" + hit.getAccession() + "', score " +
                         String(hit.getScore());
  }

} // namespace OpenMS

