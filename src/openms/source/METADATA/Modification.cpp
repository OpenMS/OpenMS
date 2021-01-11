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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Modification.h>

using namespace std;

namespace OpenMS
{

  const std::string Modification::NamesOfSpecificityType[] = {"AA", "AA_AT_CTERM", "AA_AT_NTERM", "CTERM", "NTERM"};

  Modification::Modification() :
    SampleTreatment("Modification"),
    reagent_name_(""),
    mass_(0.0),
    specificity_type_(AA),
    affected_amino_acids_("")
  {
  }

  Modification::~Modification()
  {
  }

  bool Modification::operator==(const SampleTreatment & rhs) const
  {
    if (type_ != rhs.getType())
    {
      return false;
    }

    const Modification * tmp = dynamic_cast<const Modification *>(&rhs);
    return SampleTreatment::operator==(* tmp) &&
           reagent_name_ == tmp->reagent_name_ &&
           mass_ == tmp->mass_ &&
           specificity_type_ == tmp->specificity_type_ &&
           affected_amino_acids_ == tmp->affected_amino_acids_;
  }

  SampleTreatment * Modification::clone() const
  {
    SampleTreatment * tmp = new Modification(*this);
    return tmp;
  }

  const String & Modification::getReagentName() const
  {
    return reagent_name_;
  }

  void Modification::setReagentName(const String & reagent_name)
  {
    reagent_name_ = reagent_name;
  }

  double Modification::getMass() const
  {
    return mass_;
  }

  void Modification::setMass(double mass)
  {
    mass_ = mass;
  }

  const Modification::SpecificityType & Modification::getSpecificityType() const
  {
    return specificity_type_;
  }

  void Modification::setSpecificityType(const Modification::SpecificityType & specificity_type)
  {
    specificity_type_ = specificity_type;
  }

  const String & Modification::getAffectedAminoAcids() const
  {
    return affected_amino_acids_;
  }

  void Modification::setAffectedAminoAcids(const String & affected_amino_acids)
  {
    affected_amino_acids_ = affected_amino_acids;
  }

}

