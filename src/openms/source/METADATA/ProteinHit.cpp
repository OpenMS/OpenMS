// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ProteinHit.h>


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
  ProteinHit & ProteinHit::operator=(const MetaInfoInterface & source)
  {
    MetaInfoInterface::operator=(source);
    return *this;
  }

  // equality operator
  bool ProteinHit::operator==(const ProteinHit & rhs) const
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
  bool ProteinHit::operator!=(const ProteinHit & rhs) const
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
  const String & ProteinHit::getSequence() const
  {
    return sequence_;
  }

  // returns the accession of the protein
  const String & ProteinHit::getAccession() const
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
  void ProteinHit::setSequence(const String & sequence)
  {
    sequence_ = sequence;
    sequence_.trim();
  }

  // sets the description of the protein
  void ProteinHit::setDescription(const String & description)
  {
    setMetaValue("Description", description);
  }

  // sets the accession of the protein
  void ProteinHit::setAccession(const String & accession)
  {
    accession_ = accession;
    accession_.trim();
  }

  // sets the coverage (in percent) of the protein hit based upon matched peptides
  void ProteinHit::setCoverage(const double coverage)
  {
    coverage_ = coverage;
  }

  const set<pair<Size, ResidueModification> >& ProteinHit::getModifications() const
  {
    return modifications_;
  }

  void ProteinHit::setModifications(std::set<std::pair<Size, ResidueModification> >& mods)
  {
    modifications_ = mods;
  }


} // namespace OpenMS

