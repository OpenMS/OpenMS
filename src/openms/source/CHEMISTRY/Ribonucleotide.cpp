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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace std;

namespace OpenMS
{
  const EmpiricalFormula Ribonucleotide::default_baseloss_ =
    EmpiricalFormula("C5H10O5");

  ostream& operator<<(ostream& os, const Ribonucleotide& ribo)
  {
    os << "Ribonucleotide '"
       << ribo.code_ << "' ("
       << ribo.name_ << ", "
       << ribo.formula_ << ")";
    return os;
  }

  Ribonucleotide::Ribonucleotide(
    const String& name, const String& code, const String& new_code,
    const String& html_code, const EmpiricalFormula& formula, char origin,
    double mono_mass, double avg_mass, enum TermSpecificityNuc term_spec,
    const EmpiricalFormula& baseloss_formula):
    name_(name), code_(code), new_code_(new_code), html_code_(html_code),
    formula_(formula), origin_(origin), mono_mass_(mono_mass),
    avg_mass_(avg_mass), term_spec_(term_spec),
    baseloss_formula_(baseloss_formula)
  {
  }

  Ribonucleotide::~Ribonucleotide()
  {
  }

  bool Ribonucleotide::operator==(const Ribonucleotide& ribonucleotide) const
  {
    return name_ == ribonucleotide.name_ &&
        code_ == ribonucleotide.code_ &&
        new_code_ == ribonucleotide.new_code_ &&
        html_code_ == ribonucleotide.html_code_ &&
        formula_ == ribonucleotide.formula_ &&
        origin_ == ribonucleotide.origin_ &&
        mono_mass_ == ribonucleotide.mono_mass_ &&
        avg_mass_ == ribonucleotide.avg_mass_ &&
        term_spec_ == ribonucleotide.term_spec_ &&
        baseloss_formula_ == ribonucleotide.baseloss_formula_;
  }

  const String Ribonucleotide::getCode() const
  {
    return code_;
  }

  void Ribonucleotide::setCode(const String& code)
  {
    code_ = code;
  }

  const String Ribonucleotide::getName() const
  {
    return name_;
  }

  void Ribonucleotide::setName(const String& name)
  {
    name_ = name;
  }

  double Ribonucleotide::getMonoMass() const
  {
      return mono_mass_;
  }

  void Ribonucleotide::setMonoMass(double mono_mass)
  {
    mono_mass_ = mono_mass;
  }

  double Ribonucleotide::getAvgMass() const
  {
    return avg_mass_;
  }

  void Ribonucleotide::setAvgMass(double avg_mass)
  {
    avg_mass_ = avg_mass;
  }

  const String Ribonucleotide::getNewCode() const
  {
    return new_code_;
  }

  void Ribonucleotide::setNewCode(const String& new_code)
  {
    new_code_ = new_code;
  }

  char Ribonucleotide::getOrigin() const
  {
    return origin_;
  }

  void Ribonucleotide::setOrigin(char origin)
  {
    origin_ = origin;
  }

  String Ribonucleotide::getHTMLCode() const
  {
    return html_code_;
  }

  void Ribonucleotide::setHTMLCode(const String& html_code)
  {
    html_code_ = html_code;
  }

  const EmpiricalFormula Ribonucleotide::getFormula() const
  {
    return formula_;
  }

  void Ribonucleotide::setFormula(const EmpiricalFormula& formula)
  {
    formula_ = formula;
  }

  enum Ribonucleotide::TermSpecificityNuc Ribonucleotide::getTermSpecificity() const
  {
    return term_spec_;
  }

  void Ribonucleotide::setTermSpecificity(enum TermSpecificityNuc term_spec)
  {
    if (term_spec == NUMBER_OF_TERM_SPECIFICITY)
    {
      String msg = "invalid terminal specificity";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    msg, "NUMBER_OF_TERM_SPECIFICITY");
    }
    term_spec_ = term_spec;
  }

  const EmpiricalFormula Ribonucleotide::getBaselossFormula() const
  {
    return baseloss_formula_;
  }

  void Ribonucleotide::setBaselossFormula(const EmpiricalFormula& formula)
  {
    baseloss_formula_ = formula;
  }

  bool Ribonucleotide::isModified() const
  {
    return (code_.length() != 1) || (code_[0] != origin_);
  }

  bool Ribonucleotide::isAmbiguous() const
  {
    return code_.back() == '?';
  }

}
