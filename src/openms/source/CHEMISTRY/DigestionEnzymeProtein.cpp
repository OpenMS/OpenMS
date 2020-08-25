// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
  DigestionEnzymeProtein::DigestionEnzymeProtein() :
    DigestionEnzyme(),
    n_term_gain_(""),
    c_term_gain_(""),
    psi_id_(""),
    xtandem_id_(""),
    comet_id_(-1),
    crux_id_(""),
    msgf_id_(-1),
    omssa_id_(-1)
  {
  }

  DigestionEnzymeProtein::DigestionEnzymeProtein(const DigestionEnzyme& d) :
      DigestionEnzyme(d),
      n_term_gain_(""),
      c_term_gain_(""),
      psi_id_(""),
      xtandem_id_(""),
      comet_id_(-1),
      crux_id_(""),
      msgf_id_(-1),
      omssa_id_(-1)
  {
  }

  DigestionEnzymeProtein::DigestionEnzymeProtein(const String& name,
                                                 const String& cleavage_regex,
                                                 const std::set<String>& synonyms,
                                                 String regex_description,
                                                 EmpiricalFormula n_term_gain,
                                                 EmpiricalFormula c_term_gain,
                                                 String psi_id,
                                                 String xtandem_id,
                                                 Int comet_id,
                                                 String crux_id,
                                                 Int msgf_id,
                                                 Int omssa_id) :
    DigestionEnzyme(name, cleavage_regex, synonyms, std::move(regex_description)),
    n_term_gain_(std::move(n_term_gain)),
    c_term_gain_(std::move(c_term_gain)),
    psi_id_(std::move(psi_id)),
    xtandem_id_(std::move(xtandem_id)),
    comet_id_(comet_id),
    crux_id_(std::move(crux_id)),
    msgf_id_(msgf_id),
    omssa_id_(omssa_id)
  {
  }

  DigestionEnzymeProtein::~DigestionEnzymeProtein()
  {
  }

  void DigestionEnzymeProtein::setNTermGain(const EmpiricalFormula& value)
  {
    n_term_gain_ = value;
  }

  void DigestionEnzymeProtein::setCTermGain(const EmpiricalFormula& value)
  {
    c_term_gain_ = value;
  }

  EmpiricalFormula DigestionEnzymeProtein::getNTermGain() const
  {
    return n_term_gain_;
  }

  EmpiricalFormula DigestionEnzymeProtein::getCTermGain() const
  {
    return c_term_gain_;
  }

  void DigestionEnzymeProtein::setPSIID(const String& value)
  {
    psi_id_ = value;
  }

  String DigestionEnzymeProtein::getPSIID() const
  {
    return psi_id_;
  }

  void DigestionEnzymeProtein::setXTandemID(const String& value)
  {
    xtandem_id_ = value;
  }

  String DigestionEnzymeProtein::getXTandemID() const
  {
    return xtandem_id_;
  }

  void DigestionEnzymeProtein::setCometID(Int value)
  {
    comet_id_ = value;
  }

  Int DigestionEnzymeProtein::getCometID() const
  {
    return comet_id_;
  }

  void DigestionEnzymeProtein::setCruxID(const String& value)
  {
    crux_id_ = value;
  }

  String DigestionEnzymeProtein::getCruxID() const
  {
    return crux_id_;
  }

  void DigestionEnzymeProtein::setOMSSAID(Int value)
  {
    omssa_id_ = value;
  }

  Int DigestionEnzymeProtein::getOMSSAID() const
  {
    return omssa_id_;
  }

  void DigestionEnzymeProtein::setMSGFID(Int value)
  {
    msgf_id_ = value;
  }

  Int DigestionEnzymeProtein::getMSGFID() const
  {
    return msgf_id_;
  }

  bool DigestionEnzymeProtein::operator==(const DigestionEnzymeProtein& enzyme) const
  {
    return DigestionEnzyme::operator==(enzyme) &&
           n_term_gain_ == enzyme.n_term_gain_ &&
           c_term_gain_ == enzyme.c_term_gain_ &&
           psi_id_ == enzyme.psi_id_ &&
           xtandem_id_ == enzyme.xtandem_id_ &&
           comet_id_ == enzyme.comet_id_ &&
           crux_id_ == enzyme.crux_id_ &&
           msgf_id_ == enzyme.msgf_id_ &&
           omssa_id_ == enzyme.omssa_id_;
  }

  // Note: comparison operators are not inherited. TODO rename it and make virtual
  bool DigestionEnzymeProtein::operator==(const String& cleavage_regex) const
  {
    return cleavage_regex_ == cleavage_regex;
  }

  bool DigestionEnzymeProtein::operator!=(const String& cleavage_regex) const
  {
    return cleavage_regex_ != cleavage_regex;
  }

  bool DigestionEnzymeProtein::operator!=(const DigestionEnzymeProtein& enzyme) const
  {
    return !(*this == enzyme);
  }

  bool DigestionEnzymeProtein::operator<(const DigestionEnzymeProtein& enzyme) const
  {
    return this->getName() < enzyme.getName();
  }

  bool DigestionEnzymeProtein::setValueFromFile(const String& key, const String& value)
  {
    if (DigestionEnzyme::setValueFromFile(key, value))
    {
      return true;
    }
    if (key.hasSuffix(":NTermGain"))
    {
      setNTermGain(EmpiricalFormula(value));
      return true;
    }
    if (key.hasSuffix(":CTermGain"))
    {
      setCTermGain(EmpiricalFormula(value));
      return true;
    }
    if (key.hasSuffix(":PSIID"))
    {
      setPSIID(value);
      return true;
    }
    if (key.hasSuffix(":XTandemID"))
    {
      setXTandemID(value);
      return true;
    }
    if (key.hasSuffix(":CometID"))
    {
      setCometID(value.toInt());
      return true;
    }
    if (key.hasSuffix(":CruxID"))
    {
      setCruxID(value);
      return true;
    }
    if (key.hasSuffix(":OMSSAID"))
    {
      setOMSSAID(value.toInt());
      return true;
    }
    if (key.hasSuffix(":MSGFID"))
    {
      setMSGFID(value.toInt());
      return true;
    }
    return false;
  }

  ostream& operator<<(ostream& os, const DigestionEnzymeProtein& enzyme)
  {
    os << static_cast<const DigestionEnzyme&>(enzyme) << " "
       << enzyme.psi_id_;
    return os;
  }

}

