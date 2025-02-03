// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

  Modification::~Modification() = default;

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

