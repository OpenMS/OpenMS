// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
