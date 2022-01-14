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

#include <boost/variant.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    typedef boost::variant<IdentifiedPeptideRef, IdentifiedCompoundRef,
                           IdentifiedOligoRef> RefVariant;

    struct IdentifiedMolecule: public RefVariant
    {
      IdentifiedMolecule() = default;

      IdentifiedMolecule(IdentifiedPeptideRef ref): RefVariant(ref) {};
      IdentifiedMolecule(IdentifiedCompoundRef ref): RefVariant(ref) {};
      IdentifiedMolecule(IdentifiedOligoRef ref): RefVariant(ref) {};

      IdentifiedMolecule(const IdentifiedMolecule&) = default;

      bool operator==(const IdentifiedMolecule& other) const
      {
        return RefVariant::operator==(static_cast<RefVariant>(other));
      }

      bool operator!=(const IdentifiedMolecule& other) const
      {
        return !operator==(other);
      }

      bool operator<(const IdentifiedMolecule& other) const
      {
        return RefVariant::operator<(static_cast<RefVariant>(other));
      }

      MoleculeType getMoleculeType() const
      {
        if (boost::get<IdentifiedPeptideRef>(this))
        {
          return MoleculeType::PROTEIN;
        }
        if (boost::get<IdentifiedCompoundRef>(this))
        {
          return MoleculeType::COMPOUND;
        }
        // if (boost::get<IdentifiedOligoRef>(this))
        return MoleculeType::RNA;
      }

      IdentifiedPeptideRef getIdentifiedPeptideRef() const
      {
        if (const IdentifiedPeptideRef* ref_ptr =
            boost::get<IdentifiedPeptideRef>(this))
        {
          return *ref_ptr;
        }
        String msg = "matched molecule is not a peptide";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }

      IdentifiedCompoundRef getIdentifiedCompoundRef() const
      {
        if (const IdentifiedCompoundRef* ref_ptr =
            boost::get<IdentifiedCompoundRef>(this))
        {
          return *ref_ptr;
        }
        String msg = "matched molecule is not a compound";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }

      IdentifiedOligoRef getIdentifiedOligoRef() const
      {
        if (const IdentifiedOligoRef* ref_ptr =
            boost::get<IdentifiedOligoRef>(this))
        {
          return *ref_ptr;
        }
        String msg = "matched molecule is not an oligonucleotide";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }

      String toString() const
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

      EmpiricalFormula getFormula(Size fragment_type = 0, Int charge = 0) const
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
    };
  }
}
