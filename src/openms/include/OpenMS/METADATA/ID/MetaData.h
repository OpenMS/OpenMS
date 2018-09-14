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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_ID_METADATA_H
#define OPENMS_METADATA_ID_METADATA_H

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/Software.h>

#include <boost/optional.hpp>
#include <boost/unordered_map.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /// Wrapper that adds @p operator< to iterators, so they can be used as (part of) keys in maps/sets or @p multi_index_containers
    template <typename Iterator>
    struct IteratorWrapper: public Iterator
    {
      IteratorWrapper(): Iterator() {}

      IteratorWrapper(const Iterator& it): Iterator(it) {}

      bool operator<(const IteratorWrapper& other) const
      {
        // compare by address of referenced element:
        return &(**this) < &(*other);
      }

      /// Conversion to pointer type for hashing
      operator uintptr_t() const
      {
        return uintptr_t(&(**this));
      }
    };


    enum MoleculeType
    {
      PROTEIN,
      COMPOUND,
      RNA,
      SIZE_OF_MOLECULETYPE
    };


    // Input files that were processed:
    typedef std::set<String> InputFiles;
    typedef IteratorWrapper<InputFiles::iterator> InputFileRef;


    /*!
      Information about software used for data processing.

      If the same processing is applied to multiple ID runs, e.g. if multiple files (fractions, replicates) are searched with the same search engine, store the
 software information only once.
    */
    typedef std::set<Software> DataProcessingSoftware;
    typedef IteratorWrapper<DataProcessingSoftware::iterator> ProcessingSoftwareRef;


    /*!
      Data processing step that is applied to the data (e.g. database search, PEP calculation, filtering, ConsensusID).
    */
    struct DataProcessingStep: public MetaInfoInterface
    {
      ProcessingSoftwareRef software_ref;

      std::vector<InputFileRef> input_file_refs;

      std::vector<String> primary_files; // path(s) to primary MS data

      DateTime date_time;

      // @TODO: add processing actions that are relevant for ID data
      std::set<DataProcessing::ProcessingAction> actions;

      explicit DataProcessingStep(
        ProcessingSoftwareRef software_ref,
        const std::vector<InputFileRef>& input_file_refs =
        std::vector<InputFileRef>(), const std::vector<String>& primary_files =
        std::vector<String>(), const DateTime& date_time = DateTime::now(),
        std::set<DataProcessing::ProcessingAction> actions =
        std::set<DataProcessing::ProcessingAction>()):
        software_ref(software_ref), input_file_refs(input_file_refs),
        primary_files(primary_files), date_time(date_time), actions(actions)
      {
      }

      DataProcessingStep(const DataProcessingStep& other) = default;

      // don't compare meta data (?):
      bool operator<(const DataProcessingStep& other) const
      {
        return (std::tie(software_ref, input_file_refs, primary_files,
                         date_time, actions) <
                std::tie(other.software_ref, other.input_file_refs,
                         other.primary_files, other.date_time, other.actions));
      }

      // don't compare meta data (?):
      bool operator==(const DataProcessingStep& other) const
      {
        return (std::tie(software_ref, input_file_refs, primary_files,
                         date_time, actions) ==
                std::tie(other.software_ref, other.input_file_refs,
                         other.primary_files, other.date_time, other.actions));
      }
    };

    typedef std::set<DataProcessingStep> DataProcessingSteps;
    typedef IteratorWrapper<DataProcessingSteps::iterator> ProcessingStepRef;


    enum MassType
    {
      MONOISOTOPIC,
      AVERAGE,
      SIZE_OF_MASSTYPE
    };


    /*!
      Parameters specific to a database search step.
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

    /*!
      Information about a score type.
    */
    struct ScoreType: public MetaInfoInterface
    {
      CVTerm cv_term;

      String name;

      bool higher_better;

      // reference to the software that assigned the score:
      boost::optional<ProcessingSoftwareRef> software_opt;
      // @TODO: scores assigned by different software tools/versions are
      // considered as different scores (even if they have the same name) -
      // does that make sense?

      ScoreType():
        higher_better(true), software_opt()
      {
      }

      explicit ScoreType(const CVTerm& cv_term, bool higher_better,
                         boost::optional<ProcessingSoftwareRef> software_opt =
                         boost::none):
        cv_term(cv_term), name(cv_term.getName()), higher_better(higher_better),
        software_opt(software_opt)
      {
      }

      explicit ScoreType(const String& name, bool higher_better,
                         boost::optional<ProcessingSoftwareRef> software_opt =
                         boost::none):
        cv_term(), name(name), higher_better(higher_better),
        software_opt(software_opt)
      {
      }

      ScoreType(const ScoreType& other) = default;

      // don't include "higher_better" in the comparison:
      bool operator<(const ScoreType& other) const
      {
        return (std::tie(cv_term.getAccession(), name, software_opt) <
                std::tie(other.cv_term.getAccession(), other.name,
                         other.software_opt));
      }

      // don't include "higher_better" in the comparison:
      bool operator==(const ScoreType& other) const
      {
        return (std::tie(cv_term.getAccession(), name, software_opt) ==
                std::tie(other.cv_term.getAccession(), other.name,
                         other.software_opt));
      }
    };

    typedef std::set<ScoreType> ScoreTypes;
    typedef IteratorWrapper<ScoreTypes::iterator> ScoreTypeRef;

    // @TODO: use a "boost::multi_index_container" to allow efficient access in
    // sequence and by key?
    typedef std::vector<std::pair<ScoreTypeRef, double>> ScoreList;
  }
}

#endif
