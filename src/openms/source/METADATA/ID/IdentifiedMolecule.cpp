// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ID/IdentifiedMolecule.h>

namespace OpenMS::IdentificationDataInternal
{

bool operator==(const IdentifiedMolecule& a, const IdentifiedMolecule& b)
{
  return operator==(static_cast<RefVariant>(a), static_cast<RefVariant>(b));
}

bool operator!=(const IdentifiedMolecule& a, const IdentifiedMolecule& b)
{
  return !operator==(a, b);
}

bool operator<(const IdentifiedMolecule& a, const IdentifiedMolecule& b)
{
  return operator<(static_cast<RefVariant>(a), static_cast<RefVariant>(b));
}

MoleculeType IdentifiedMolecule::getMoleculeType() const
{
  if (std::get_if<IdentifiedPeptideRef>(this))
  {
    return MoleculeType::PROTEIN;
  }
  if (std::get_if<IdentifiedCompoundRef>(this))
  {
    return MoleculeType::COMPOUND;
  }
  // if (get<IdentifiedOligoRef>(this))
  return MoleculeType::RNA;
}

IdentifiedPeptideRef IdentifiedMolecule::getIdentifiedPeptideRef() const
{
  if (const IdentifiedPeptideRef* ref_ptr =
      std::get_if<IdentifiedPeptideRef>(this))
  {
    return *ref_ptr;
  }
  String msg = "matched molecule is not a peptide";
  throw Exception::IllegalArgument(__FILE__, __LINE__,
                                    OPENMS_PRETTY_FUNCTION, msg);
}

IdentifiedCompoundRef IdentifiedMolecule::getIdentifiedCompoundRef() const
{
  if (const IdentifiedCompoundRef* ref_ptr =
      std::get_if<IdentifiedCompoundRef>(this))
  {
    return *ref_ptr;
  }
  String msg = "matched molecule is not a compound";
  throw Exception::IllegalArgument(__FILE__, __LINE__,
                                    OPENMS_PRETTY_FUNCTION, msg);
}


IdentifiedOligoRef IdentifiedMolecule::getIdentifiedOligoRef() const
{
  if (const IdentifiedOligoRef* ref_ptr =
      std::get_if<IdentifiedOligoRef>(this))
  {
    return *ref_ptr;
  }
  String msg = "matched molecule is not an oligonucleotide";
  throw Exception::IllegalArgument(__FILE__, __LINE__,
                                    OPENMS_PRETTY_FUNCTION, msg);
}

EmpiricalFormula IdentifiedMolecule::getFormula(Size fragment_type /*= 0*/, Int charge /*= 0*/) const
{
  switch (getMoleculeType())
  {
    case MoleculeType::PROTEIN:
    {
      auto type = static_cast<Residue::ResidueType>(fragment_type);
      return getIdentifiedPeptideRef()->sequence.getFormula(type, charge);
    }
    case MoleculeType::COMPOUND:
    {
      // @TODO: what about fragment type and charge?
      return getIdentifiedCompoundRef()->formula;
    }
    case MoleculeType::RNA:
    {
      auto type = static_cast<NASequence::NASFragmentType>(fragment_type);
      return getIdentifiedOligoRef()->sequence.getFormula(type, charge);
    }
    default:
      throw Exception::NotImplemented(__FILE__, __LINE__,
                                      OPENMS_PRETTY_FUNCTION);
  }
}

String IdentifiedMolecule::toString() const
{
  switch (getMoleculeType())
  {
    case MoleculeType::PROTEIN:
      return getIdentifiedPeptideRef()->sequence.toString();
    case MoleculeType::COMPOUND:
      return getIdentifiedCompoundRef()->identifier; // or use "name"?
    case MoleculeType::RNA:
      return getIdentifiedOligoRef()->sequence.toString();
    default:
      throw Exception::NotImplemented(__FILE__, __LINE__,
                                      OPENMS_PRETTY_FUNCTION);
  }
}

} // namespace
