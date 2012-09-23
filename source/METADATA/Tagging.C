// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Tagging.h>

using namespace std;

namespace OpenMS
{
  const std::string Tagging::NamesOfIsotopeVariant[] = {"LIGHT", "HEAVY"};

  Tagging::Tagging() :
    Modification(),
    mass_shift_(0.0),
    variant_(LIGHT)
  {
    type_ = "Tagging";
  }

  Tagging::Tagging(const Tagging & source) :
    Modification(source),
    mass_shift_(source.mass_shift_),
    variant_(source.variant_)
  {
    //????
  }

  Tagging::~Tagging()
  {
    //????
  }

  Tagging & Tagging::operator=(const Tagging & source)
  {
    if (&source == this)
      return *this;

    Modification::operator=(source);
    mass_shift_ = source.mass_shift_;
    variant_ = source.variant_;

    return *this;
  }

  bool Tagging::operator==(const SampleTreatment & rhs) const
  {
    if (type_ != rhs.getType())
      return false;

    const Tagging * tmp = dynamic_cast<const Tagging *>(&rhs);
    return Modification::operator==(rhs) &&
           mass_shift_ == tmp->mass_shift_ &&
           variant_ == tmp->variant_
    ;
  }

  SampleTreatment * Tagging::clone() const
  {
    SampleTreatment * tmp = new Tagging(*this);
    return tmp;
  }

  DoubleReal Tagging::getMassShift() const
  {
    return mass_shift_;
  }

  void Tagging::setMassShift(DoubleReal mass_shift)
  {
    mass_shift_ = mass_shift;
  }

  const Tagging::IsotopeVariant & Tagging::getVariant() const
  {
    return variant_;
  }

  void Tagging::setVariant(const Tagging::IsotopeVariant & variant)
  {
    variant_ = variant;
  }

}
