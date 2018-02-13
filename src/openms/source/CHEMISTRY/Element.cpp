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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/Element.h>

#include <ostream>

using namespace std;

namespace OpenMS
{
  Element::Element() :
    name_(OPENMS_CHEMISTRY_ELEMENT_NAME_DEFAULT),
    symbol_(OPENMS_CHEMISTRY_ELEMENT_SYMBOL_DEFAULT),
    atomic_number_(OPENMS_CHEMISTRY_ELEMENT_ATOMICNUMBER_DEFAULT),
    average_weight_(OPENMS_CHEMISTRY_ELEMENT_WEIGHT_DEFAULT),
    mono_weight_(OPENMS_CHEMISTRY_ELEMENT_WEIGHT_DEFAULT)
  {
  }

  Element::Element(const Element & e) :
    name_(e.name_),
    symbol_(e.symbol_),
    atomic_number_(e.atomic_number_),
    average_weight_(e.average_weight_),
    mono_weight_(e.mono_weight_),
    isotopes_(e.isotopes_)
  {
  }

  Element::Element(const String & name,
                   const String & symbol,
                   UInt atomic_number,
                   double average_weight,
                   double mono_weight,
                   const IsotopeDistribution & isotopes) :
    name_(name),
    symbol_(symbol),
    atomic_number_(atomic_number),
    average_weight_(average_weight),
    mono_weight_(mono_weight),
    isotopes_(isotopes)
  {
  }

  Element::~Element()
  {
  }

  void Element::setAtomicNumber(UInt atomic_number)
  {
    atomic_number_ = atomic_number;
  }

  UInt Element::getAtomicNumber() const
  {
    return atomic_number_;
  }

  void Element::setAverageWeight(double weight)
  {
    average_weight_ = weight;
  }

  double Element::getAverageWeight() const
  {
    return average_weight_;
  }

  void Element::setMonoWeight(double weight)
  {
    mono_weight_ = weight;
  }

  double Element::getMonoWeight() const
  {
    return mono_weight_;
  }

  void Element::setIsotopeDistribution(const IsotopeDistribution & distribution)
  {
    isotopes_ = distribution;
  }

  const IsotopeDistribution & Element::getIsotopeDistribution() const
  {
    return isotopes_;
  }

  void Element::setName(const String & name)
  {
    name_ = name;
  }

  const String & Element::getName() const
  {
    return name_;
  }

  void Element::setSymbol(const String & symbol)
  {
    symbol_ = symbol;
  }

  const String & Element::getSymbol() const
  {
    return symbol_;
  }

  Element & Element::operator=(const Element & element)
  {
    name_ = element.name_;
    symbol_ = element.symbol_;
    atomic_number_ = element.atomic_number_;
    average_weight_ = element.average_weight_;
    mono_weight_ = element.mono_weight_;
    isotopes_ = element.isotopes_;
    return *this;
  }

  bool Element::operator==(const Element & element) const
  {
    return name_ == element.name_ &&
           symbol_ == element.symbol_ &&
           atomic_number_ == element.atomic_number_ &&
           average_weight_ == element.average_weight_ &&
           mono_weight_ == element.mono_weight_ &&
           isotopes_ == element.isotopes_;
  }

  bool Element::operator!=(const Element & element) const
  {
    return !(*this == element);
  }

  std::ostream & operator<<(std::ostream & os, const Element & element)
  {
    os  << element.name_ << " "
    << element.symbol_ << " "
    << element.atomic_number_ << " "
    << element.average_weight_ << " "
    << element.mono_weight_;

    for (auto const & isotope : element.isotopes_)
    {
      if (isotope.getIntensity() > 0.0f)
      {
        os << " " << isotope.getPosition() << "=" << isotope.getIntensity() * 100 << "%";
      }
    }
    return os;
  }

} // namespace OpenMS
