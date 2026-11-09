// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
  /**
    @brief Settings of a database search: the database, charges, modifications, tolerances and digestion

    Shared by the established classes (ProteinIdentification::SearchParameters is this class) and the
    run settings of IdentificationData.

    @ingroup Metadata
  */
  struct OPENMS_DLLAPI SearchParameters :
    public MetaInfoInterface
  {
    /// Peak mass type
    enum class PeakMassType
    {
      MONOISOTOPIC,
      AVERAGE,
      SIZE_OF_PEAKMASSTYPE
    };

    std::string db; ///< The used database
    std::string db_version; ///< The database version
    std::string taxonomy; ///< The taxonomy restriction
    std::string charges; ///< The allowed charges for the search
    PeakMassType mass_type; ///< Mass type of the peaks
    std::vector<std::string> fixed_modifications; ///< Used fixed modifications
    std::vector<std::string> variable_modifications; ///< Allowed variable modifications
    UInt missed_cleavages; ///< The number of allowed missed cleavages
    double fragment_mass_tolerance; ///< Mass tolerance of fragment ions (Dalton or ppm)
    bool fragment_mass_tolerance_ppm; ///< Mass tolerance unit of fragment ions (true: ppm, false: Dalton)
    double precursor_mass_tolerance; ///< Mass tolerance of precursor ions (Dalton or ppm)
    bool precursor_mass_tolerance_ppm; ///< Mass tolerance unit of precursor ions (true: ppm, false: Dalton)
    Protease digestion_enzyme; ///< The cleavage site information in details (from ProteaseDB)
    EnzymaticDigestion::Specificity enzyme_term_specificity; ///< The number of required cutting-rule matching termini during search (none=0, semi=1, or full=2)

    SearchParameters();
    /// Copy constructor
    SearchParameters(const SearchParameters&) = default;
    /// Move constructor
    SearchParameters(SearchParameters&&) = default;
    /// Destructor
    ~SearchParameters() = default;

    /// Assignment operator
    SearchParameters& operator=(const SearchParameters&) = default;
    /// Move assignment operator
    SearchParameters& operator=(SearchParameters&&)& = default;

    bool operator==(const SearchParameters& rhs) const;

    bool operator!=(const SearchParameters& rhs) const;

    /// returns the charge range from the search engine settings as a pair of ints
    std::pair<int,int> getChargeRange() const;

    /// Tests if these search engine settings are mergeable with @p sp
    /// depending on the given @p experiment_type.
    /// Modifications are compared as sets. Databases based on filename.
    /// "labeled_MS1" experiments additionally allow different modifications.
    bool mergeable(const SearchParameters& sp, const std::string& experiment_type) const;

    private:
    int getChargeValue_(std::string& charge_str) const;
  };
} // namespace OpenMS
