// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Xiao Liang  $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <cstdlib>
#include <iostream>

using namespace std;

namespace OpenMS
{
  DigestionEnzyme::DigestionEnzyme() :
    name_("unknown_enzyme"),
    cleavage_regex_(""),
    synonyms_(),
    regex_description_(""),
    n_term_gain_(""),
    c_term_gain_(""),
    psi_id_(""),
    xtandem_id_(""),
    comet_id_(),
    msgf_id_(-1),
    omssa_id_()
  {
  }

  DigestionEnzyme::DigestionEnzyme(const DigestionEnzyme& enzyme) :
    name_(enzyme.name_),
    cleavage_regex_(enzyme.cleavage_regex_),
    synonyms_(enzyme.synonyms_),
    regex_description_(enzyme.regex_description_),
    n_term_gain_(enzyme.n_term_gain_),
    c_term_gain_(enzyme.c_term_gain_),
    psi_id_(enzyme.psi_id_),
    xtandem_id_(enzyme.xtandem_id_),
    comet_id_(enzyme.comet_id_),
    msgf_id_(enzyme.msgf_id_),
    omssa_id_(enzyme.omssa_id_)
  {
  }

  DigestionEnzyme::DigestionEnzyme(const String& name,
                                   const String& cleavage_regex,
                                   const std::set<String>& synonyms,
                                   String regex_description,
                                   EmpiricalFormula n_term_gain,
                                   EmpiricalFormula c_term_gain,
                                   String psi_id,
                                   String xtandem_id,
                                   UInt comet_id,
                                   Int msgf_id,
                                   UInt omssa_id) :
    name_(name),
    cleavage_regex_(cleavage_regex),
    synonyms_(synonyms),
    regex_description_(regex_description),
    n_term_gain_(n_term_gain),
    c_term_gain_(c_term_gain),
    psi_id_(psi_id),
    xtandem_id_(xtandem_id),
    comet_id_(comet_id),
    msgf_id_(msgf_id),
    omssa_id_(omssa_id)
  {
  }

  DigestionEnzyme::~DigestionEnzyme()
  {
  }

  DigestionEnzyme& DigestionEnzyme::operator=(const DigestionEnzyme& enzyme)
  {
    if (this != &enzyme)
    {
      name_ = enzyme.name_;
      cleavage_regex_ = enzyme.cleavage_regex_;
      synonyms_ = enzyme.synonyms_;
      regex_description_ = enzyme.regex_description_;
      n_term_gain_ = enzyme.n_term_gain_;
      c_term_gain_ = enzyme.c_term_gain_;
      psi_id_ = enzyme.psi_id_;
      xtandem_id_ = enzyme.xtandem_id_;
      comet_id_ = enzyme.comet_id_;
      omssa_id_ = enzyme.omssa_id_;
      msgf_id_ = enzyme.msgf_id_;
    }
    return *this;
  }

  void DigestionEnzyme::setName(const String& name)
  {
    name_ = name;
  }

  const String& DigestionEnzyme::getName() const
  {
    return name_;
  }

  void DigestionEnzyme::setSynonyms(const set<String>& synonyms)
  {
    synonyms_ = synonyms;
  }

  void DigestionEnzyme::addSynonym(const String& synonym)
  {
    synonyms_.insert(synonym);
  }

  const set<String>& DigestionEnzyme::getSynonyms() const
  {
    return synonyms_;
  }

  void DigestionEnzyme::setRegEx(const String& cleavage_regex)
  {
    cleavage_regex_ = cleavage_regex;
  }

  const String& DigestionEnzyme::getRegEx() const
  {
    return cleavage_regex_;
  }

  void DigestionEnzyme::setRegExDescription(String value)
  {
    regex_description_ = value;
  }

  String DigestionEnzyme::getRegExDescription() const
  {
    return regex_description_;
  }

  void DigestionEnzyme::setNTermGain(EmpiricalFormula value)
  {
    n_term_gain_ = value;
  }

  void DigestionEnzyme::setCTermGain(EmpiricalFormula value)
  {
    c_term_gain_ = value;
  }

  EmpiricalFormula DigestionEnzyme::getNTermGain() const
  {
    return n_term_gain_;
  }

  EmpiricalFormula DigestionEnzyme::getCTermGain() const
  {
    return c_term_gain_;
  }

  void DigestionEnzyme::setPSIID(String value)
  {
    psi_id_ = value;
  }

  String DigestionEnzyme::getPSIID() const
  {
    return psi_id_;
  }

  void DigestionEnzyme::setXTandemID(String value)
  {
    xtandem_id_ = value;
  }

  String DigestionEnzyme::getXTandemID() const
  {
    return xtandem_id_;
  }

  void DigestionEnzyme::setCometID(UInt value)
  {
    comet_id_ = value;
  }

  UInt DigestionEnzyme::getCometID() const
  {
    return comet_id_;
  }

  void DigestionEnzyme::setOMSSAID(UInt value)
  {
    omssa_id_ = value;
  }

  UInt DigestionEnzyme::getOMSSAID() const
  {
    return omssa_id_;
  }

  void DigestionEnzyme::setMSGFID(Int value)
  {
    msgf_id_ = value;
  }

  Int DigestionEnzyme::getMSGFID() const
  {
    return msgf_id_;
  }

  bool DigestionEnzyme::operator==(const DigestionEnzyme& enzyme) const
  {
    return name_ == enzyme.name_ &&
           synonyms_ == enzyme.synonyms_ &&
           cleavage_regex_ == enzyme.cleavage_regex_ &&
           regex_description_ == enzyme.regex_description_ &&
           n_term_gain_ == enzyme.n_term_gain_ &&
           c_term_gain_ == enzyme.c_term_gain_ &&
           psi_id_ == enzyme.psi_id_ &&
           xtandem_id_ == enzyme.xtandem_id_ &&
           comet_id_ == enzyme.comet_id_ &&
           msgf_id_ == enzyme.msgf_id_ &&
           omssa_id_ == enzyme.omssa_id_;
  }

  bool DigestionEnzyme::operator==(String cleavage_regex) const
  {
    return cleavage_regex_ == cleavage_regex;
  }

  bool DigestionEnzyme::operator!=(String cleavage_regex) const
  {
    return cleavage_regex_ != cleavage_regex;
  }

  bool DigestionEnzyme::operator!=(const DigestionEnzyme& enzyme) const
  {
    return !(*this == enzyme);
  }

  bool DigestionEnzyme::operator<(const DigestionEnzyme& enzyme) const
  {
    return this->getName() < enzyme.getName();
  }

  ostream& operator<<(ostream& os, const DigestionEnzyme& enzyme)
  {
    os << enzyme.name_ << " "
    << enzyme.cleavage_regex_ << " "
    << enzyme.regex_description_ << " "
    << enzyme.psi_id_;
    return os;
  }

}
