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

#include <OpenMS/CHEMISTRY/Enzyme.h>
#include <cstdlib>
#include <iostream>

using namespace std;

namespace OpenMS
{
  // enzyme
  Enzyme::Enzyme() :
    name_("unknown_enzyme"),
    cleavage_regex_(""),
    synonyms_(),
    regex_description_(""),
    n_term_gain_(""),
    c_term_gain_(""),
    psi_id_(""),
    xtandem_id_(""),
    omssa_id_()
  {
  }

  Enzyme::Enzyme(const Enzyme & enzyme) :
    name_(enzyme.name_),
    cleavage_regex_(enzyme.cleavage_regex_),
    synonyms_(enzyme.synonyms_),
    regex_description_(enzyme.regex_description_),
    n_term_gain_(enzyme.n_term_gain_),
    c_term_gain_(enzyme.c_term_gain_),
    psi_id_(enzyme.psi_id_),
    xtandem_id_(enzyme.xtandem_id_),
    omssa_id_(enzyme.omssa_id_)
  {
  }

  Enzyme::Enzyme(const String & name,
                   const String & cleavage_regex,
                   const std::set<String> & synonyms,
                   String regex_description,
                   EmpiricalFormula n_term_gain,
                   EmpiricalFormula c_term_gain,
                   String psi_id,
                   String xtandem_id,
                   UInt omssa_id) :
    name_(name),
    cleavage_regex_(cleavage_regex),
    synonyms_(synonyms),
    regex_description_(regex_description),
    n_term_gain_(n_term_gain),
    c_term_gain_(c_term_gain),
    psi_id_(psi_id),
    xtandem_id_(xtandem_id),
    omssa_id_(omssa_id)
  {
  }
  
  Enzyme::~Enzyme()
  {
  }

  Enzyme & Enzyme::operator=(const Enzyme & enzyme)
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
      omssa_id_ = enzyme.omssa_id_;
    }
    return *this;
  }

  void Enzyme::setName(const String & name)
  {
    name_ = name;
  }

  const String & Enzyme::getName() const
  {
    return name_;
  }

  void Enzyme::setSynonyms(const set<String> & synonyms)
  {
    synonyms_ = synonyms;
  }

  void Enzyme::addSynonym(const String & synonym)
  {
    synonyms_.insert(synonym);
  }

  const set<String> & Enzyme::getSynonyms() const
  {
    return synonyms_;
  }

  void Enzyme::setRegEx(const String & cleavage_regex)
  {
    cleavage_regex_ = cleavage_regex;
  }

  const String & Enzyme::getRegEx() const
  {
    return cleavage_regex_;
  }
  
  void Enzyme::setRegExDescription(String value)
  {
    regex_description_ = value;
  }
  
  String Enzyme::getRegExDescription() const
  {
    return regex_description_;
  }

  void Enzyme::setNTermGain(EmpiricalFormula value)
  {
    n_term_gain_ = value;
  }

  void Enzyme::setCTermGain(EmpiricalFormula value)
  {
    c_term_gain_ = value;
  }
  
  EmpiricalFormula Enzyme::getNTermGain() const
  {
    return n_term_gain_;
  }

  EmpiricalFormula Enzyme::getCTermGain() const
  {
    return c_term_gain_;
  }

  void Enzyme::setPSIid(String value)
  {
    psi_id_ = value;
  }
  
  String Enzyme::getPSIid() const
  {
    return psi_id_;
  }

  void Enzyme::setXTANDEMid(String value)
  {
    xtandem_id_ = value;
  }

  String Enzyme::getXTANDEMid() const
  {
    return xtandem_id_;
  }

  void Enzyme::setOMSSAid(UInt value)
  {
    omssa_id_ = value;
  }
  
  UInt Enzyme::getOMSSAid() const
  {
    return omssa_id_;
  }
  
  bool Enzyme::operator==(const Enzyme & enzyme) const
  {
    return name_ == enzyme.name_ &&
           synonyms_ == enzyme.synonyms_ &&
           cleavage_regex_ == enzyme.cleavage_regex_ &&
           regex_description_ == enzyme.regex_description_ &&
           n_term_gain_ == enzyme.n_term_gain_ &&
           c_term_gain_ == enzyme.c_term_gain_ &&
           psi_id_ == enzyme.psi_id_ &&
           xtandem_id_ == enzyme.xtandem_id_ &&
           omssa_id_ == enzyme.omssa_id_;
  }

  bool Enzyme::operator==(String cleavage_regex) const
  {
    return cleavage_regex_ == cleavage_regex;
  }

  bool Enzyme::operator!=(String cleavage_regex) const
  {
    return cleavage_regex_ != cleavage_regex;
  }

  bool Enzyme::operator!=(const Enzyme & enzyme) const
  {
    return !(*this == enzyme);
  }

  bool Enzyme::operator<(const Enzyme & enzyme) const
  {
    return this->getName() < enzyme.getName();
  }

  ostream & operator<<(ostream & os, const Enzyme & enzyme)
  {
    os << enzyme.name_ << " "
    << enzyme.cleavage_regex_ << " "
    << enzyme.regex_description_ << " "
    << enzyme.psi_id_;
    return os;
  }

}
