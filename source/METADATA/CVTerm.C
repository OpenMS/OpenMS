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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/CVTerm.h>

using namespace std;

namespace OpenMS
{
  // CV term implementation
  CVTerm::CVTerm()
  {
  }

  CVTerm::CVTerm(const String & accession, const String & name, const String & cv_identifier_ref, const String & value, const Unit & unit) :
    accession_(accession),
    name_(name),
    cv_identifier_ref_(cv_identifier_ref),
    unit_(unit),
    value_(value)
  {
  }

  CVTerm::CVTerm(const CVTerm & rhs) :
    accession_(rhs.accession_),
    name_(rhs.name_),
    cv_identifier_ref_(rhs.cv_identifier_ref_),
    unit_(rhs.unit_),
    value_(rhs.value_)
  {
  }

  CVTerm::~CVTerm()
  {
  }

  CVTerm & CVTerm::operator=(const CVTerm & rhs)
  {
    if (this != &rhs)
    {
      accession_ = rhs.accession_;
      name_ = rhs.name_;
      cv_identifier_ref_ = rhs.cv_identifier_ref_;
      unit_ = rhs.unit_;
      value_ = rhs.value_;
    }
    return *this;
  }

  bool CVTerm::operator==(const CVTerm & rhs) const
  {
    return accession_ == rhs.accession_ &&
           name_ == rhs.name_ &&
           cv_identifier_ref_ == rhs.cv_identifier_ref_ &&
           unit_ == rhs.unit_ &&
           value_ == rhs.value_;
  }

  bool CVTerm::operator!=(const CVTerm & rhs) const
  {
    return !(*this == rhs);
  }

  void CVTerm::setAccession(const String & accession)
  {
    accession_ = accession;
  }

  const String & CVTerm::getAccession() const
  {
    return accession_;
  }

  void CVTerm::setName(const String & name)
  {
    name_ = name;
  }

  const String & CVTerm::getName() const
  {
    return name_;
  }

  void CVTerm::setCVIdentifierRef(const String & cv_identifier_ref)
  {
    cv_identifier_ref_ = cv_identifier_ref;
  }

  const String & CVTerm::getCVIdentifierRef() const
  {
    return cv_identifier_ref_;
  }

  void CVTerm::setUnit(const Unit & unit)
  {
    unit_ = unit;
  }

  const CVTerm::Unit & CVTerm::getUnit() const
  {
    return unit_;
  }

  void CVTerm::setValue(const DataValue & value)
  {
    value_ = value;
  }

  const DataValue & CVTerm::getValue() const
  {
    return value_;
  }

  bool CVTerm::hasUnit() const
  {
    return unit_.accession != "";
  }

  bool CVTerm::hasValue() const
  {
    return value_ != DataValue::EMPTY;
  }

} // namespace OpenMS
