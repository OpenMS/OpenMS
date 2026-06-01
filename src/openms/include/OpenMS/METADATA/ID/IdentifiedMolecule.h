// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/MetaData.h>
#include <OpenMS/METADATA/ID/IdentifiedCompound.h>
#include <OpenMS/METADATA/ID/IdentifiedSequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <variant>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    typedef std::variant<IdentifiedPeptideRef, IdentifiedCompoundRef,
                           IdentifiedOligoRef> RefVariant;

    /**
      @brief Variant type holding Peptide/Compound/Oligo references and convenience functions.
    **/
    struct OPENMS_DLLAPI IdentifiedMolecule: public RefVariant
    {
      IdentifiedMolecule() = default;

      IdentifiedMolecule(IdentifiedPeptideRef ref): RefVariant(ref) {};
      IdentifiedMolecule(IdentifiedCompoundRef ref): RefVariant(ref) {};
      IdentifiedMolecule(IdentifiedOligoRef ref): RefVariant(ref) {};

      IdentifiedMolecule(const IdentifiedMolecule&) = default;

      MoleculeType getMoleculeType() const;

      IdentifiedPeptideRef getIdentifiedPeptideRef() const;

      IdentifiedCompoundRef getIdentifiedCompoundRef() const;

      IdentifiedOligoRef getIdentifiedOligoRef() const;

      String toString() const;

      EmpiricalFormula getFormula(Size fragment_type = 0, Int charge = 0) const;
    };

    OPENMS_DLLAPI bool operator==(const IdentifiedMolecule& a, const IdentifiedMolecule& b);

    OPENMS_DLLAPI bool operator!=(const IdentifiedMolecule& a, const IdentifiedMolecule& b);

    OPENMS_DLLAPI bool operator<(const IdentifiedMolecule& a, const IdentifiedMolecule& b);

  }
}
