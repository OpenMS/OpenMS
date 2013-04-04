// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideEvidence.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{
  // default constructor
  PeptideEvidence::PeptideEvidence() :
    MetaInfoInterface(),
    db_sequence_ref_(""),
    translation_table_ref_(""),
    start_(-1),
    end_(-1),
    pre_('?'),
    post_('?'),
    id_(""),
    name_(""),
    missed_cleavages_(-1),
    is_decoy_(false),
    frame_(-1)
  {
  }

  // copy constructor
  PeptideEvidence::PeptideEvidence(const PeptideEvidence & rhs) :
    MetaInfoInterface(rhs),
    db_sequence_ref_(rhs.db_sequence_ref_),
    translation_table_ref_(rhs.translation_table_ref_),
    start_(rhs.start_),
    end_(rhs.end_),
    pre_(rhs.pre_),
    post_(rhs.post_),
    id_(rhs.id_),
    name_(rhs.name_),
    missed_cleavages_(rhs.missed_cleavages_),
    is_decoy_(rhs.is_decoy_),
    frame_(rhs.frame_)
  {
  }

  // destructor
  PeptideEvidence::~PeptideEvidence()
  {
  }

  PeptideEvidence & PeptideEvidence::operator=(const PeptideEvidence & rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }

    MetaInfoInterface::operator=(rhs);
    db_sequence_ref_ = rhs.db_sequence_ref_;
    translation_table_ref_ = rhs.translation_table_ref_;
    start_ = rhs.start_;
    end_ = rhs.end_;
    pre_ = rhs.pre_;
    post_ = rhs.post_;
    id_ = rhs.id_;
    name_ = name_;
    missed_cleavages_ = rhs.missed_cleavages_;
    is_decoy_ = rhs.is_decoy_;
    frame_ = rhs.frame_;

    return *this;
  }

  bool PeptideEvidence::operator==(const PeptideEvidence & rhs) const
  {
    return MetaInfoInterface::operator==(rhs) &&
           db_sequence_ref_ == rhs.db_sequence_ref_ &&
           translation_table_ref_ == rhs.translation_table_ref_ &&
           start_ == rhs.start_ &&
           end_ == rhs.end_ &&
           pre_ == rhs.pre_ &&
           post_ == rhs.post_ &&
           id_ == rhs.id_ &&
           name_ == rhs.name_ &&
           missed_cleavages_ == rhs.missed_cleavages_ &&
           is_decoy_ == rhs.is_decoy_ &&
           frame_ == rhs.frame_;
  }

  bool PeptideEvidence::operator!=(const PeptideEvidence & rhs) const
  {
    return !operator==(rhs);
  }

  const String & PeptideEvidence::getDBSequenceRef() const
  {
    return db_sequence_ref_;
  }

  void PeptideEvidence::setDBSequenceRef(const String & rhs)
  {
    db_sequence_ref_ = rhs;
  }

  const String & PeptideEvidence::getTranslationTableRef() const
  {
    return translation_table_ref_;
  }

  void PeptideEvidence::setTranslationTableRef(const String & rhs)
  {
    translation_table_ref_ = rhs;
  }

  void PeptideEvidence::setStart(Int start)
  {
    start_ = start;
  }

  Int PeptideEvidence::getStart() const
  {
    return start_;
  }

  void PeptideEvidence::setEnd(Int end)
  {
    end_ = end;
  }

  Int PeptideEvidence::getEnd() const
  {
    return end_;
  }

  void PeptideEvidence::setPre(char rhs)
  {
    pre_ = rhs;
  }

  char PeptideEvidence::getPre() const
  {
    return pre_;
  }

  void PeptideEvidence::setPost(char rhs)
  {
    post_ = rhs;
  }

  char PeptideEvidence::getPost() const
  {
    return post_;
  }

  void PeptideEvidence::setId(const String & id)
  {
    id_ = id;
  }

  const String & PeptideEvidence::getId() const
  {
    return id_;
  }

  void PeptideEvidence::setName(const String & name)
  {
    name_ = name;
  }

  const String & PeptideEvidence::getName() const
  {
    return name_;
  }

  void PeptideEvidence::setMissedCleavages(Int rhs)
  {
    missed_cleavages_ = rhs;
  }

  Int PeptideEvidence::getMissedCleavages() const
  {
    return missed_cleavages_;
  }

  void PeptideEvidence::setIsDecoy(bool is_decoy)
  {
    is_decoy_ = is_decoy;
  }

  bool PeptideEvidence::getIsDecoy() const
  {
    return is_decoy_;
  }

  void PeptideEvidence::setFrame(Int frame)
  {
    frame_ = frame;
  }

  Int PeptideEvidence::getFrame() const
  {
    return frame_;
  }

} // namespace OpenMS
