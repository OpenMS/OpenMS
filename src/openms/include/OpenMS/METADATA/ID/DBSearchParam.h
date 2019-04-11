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

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/METADATA/ID/MetaData.h>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Parameters specific to a database search step.
    */
    struct DBSearchParam: public MetaInfoInterface
    {
      enum MoleculeType molecule_type;
      enum MassType mass_type;

      String database;
      String database_version;
      String taxonomy;

      std::set<Int> charges;

      std::set<String> fixed_mods;
      std::set<String> variable_mods;

      double precursor_mass_tolerance;
      double fragment_mass_tolerance;
      bool precursor_tolerance_ppm;
      bool fragment_tolerance_ppm;

      // allow for either "DigestionEnzymeProtein" or "DigestionEnzymeRNA":
      const DigestionEnzyme* digestion_enzyme;
      Size missed_cleavages;
      Size min_length;
      Size max_length;

      DBSearchParam():
        molecule_type(MoleculeType::PROTEIN),
        mass_type(MassType::MONOISOTOPIC),
        precursor_mass_tolerance(0.0), fragment_mass_tolerance(0.0),
        precursor_tolerance_ppm(false), fragment_tolerance_ppm(false),
        digestion_enzyme(0), missed_cleavages(0), min_length(0), max_length(0)
      {
      }

      DBSearchParam(const DBSearchParam& other) = default;

      bool operator<(const DBSearchParam& other) const
      {
        return (std::tie(molecule_type, mass_type, database,
                         database_version, taxonomy, charges, fixed_mods,
                         variable_mods, fragment_mass_tolerance,
                         precursor_mass_tolerance, fragment_tolerance_ppm,
                         precursor_tolerance_ppm, digestion_enzyme,
                         missed_cleavages, min_length, max_length) <
                std::tie(other.molecule_type, other.mass_type,
                         other.database, other.database_version, other.taxonomy,
                         other.charges, other.fixed_mods, other.variable_mods,
                         other.fragment_mass_tolerance,
                         other.precursor_mass_tolerance,
                         other.fragment_tolerance_ppm,
                         other.precursor_tolerance_ppm,
                         other.digestion_enzyme, other.missed_cleavages,
                         other.min_length, other.max_length));
      }

      bool operator==(const DBSearchParam& other) const
      {
        return (std::tie(molecule_type, mass_type, database,
                         database_version, taxonomy, charges, fixed_mods,
                         variable_mods, fragment_mass_tolerance,
                         precursor_mass_tolerance, fragment_tolerance_ppm,
                         precursor_tolerance_ppm, digestion_enzyme,
                         missed_cleavages, min_length, max_length) ==
                std::tie(other.molecule_type, other.mass_type,
                         other.database, other.database_version, other.taxonomy,
                         other.charges, other.fixed_mods, other.variable_mods,
                         other.fragment_mass_tolerance,
                         other.precursor_mass_tolerance,
                         other.fragment_tolerance_ppm,
                         other.precursor_tolerance_ppm,
                         other.digestion_enzyme, other.missed_cleavages,
                         other.min_length, other.max_length));
      }
    };

    typedef std::set<DBSearchParam> DBSearchParams;
    typedef IteratorWrapper<DBSearchParams::iterator> SearchParamRef;
    typedef std::map<ProcessingStepRef, SearchParamRef> DBSearchSteps;

  }
}
