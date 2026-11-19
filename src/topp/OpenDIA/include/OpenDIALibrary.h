// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS::OpenDIA
{
  /// One library fragment (a PQP TRANSITION row with its precursor's m/z and 1/K0).
  struct LibraryTransition
  {
    std::string transition_name;   ///< PQP TRANSITION.ID.
    std::string peptide_ref;       ///< PQP PRECURSOR.ID.
    double library_intensity = 0.0;
    double product_mz = 0.0;
    double precursor_mz = 0.0;
    double precursor_im = -1.0;    ///< Library 1/K0, -1 if absent.
    std::int8_t fragment_charge = 0;
    std::int16_t fragment_nr = -1;
    std::string fragment_type;     ///< Normalised ion type ("b", "y", "b-H2O", "prec", "" none, "unknown").
    bool decoy = false;

    bool getDecoy() const { return decoy; }
    void setDecoy(bool d) { decoy = d; }
    /// "y7^2", "b3", "prec", or "" without an ion type.
    std::string getAnnotation() const;
  };

  /// One library precursor.
  struct LibraryCompound
  {
    std::string id;                ///< PQP PRECURSOR.ID.
    double rt = std::numeric_limits<double>::quiet_NaN();  ///< PRECURSOR.LIBRARY_RT (0 when NULL in the file).
    double drift_time = -1.0;      ///< PRECURSOR.LIBRARY_DRIFT_TIME, -1 if absent.
    int charge = 0;
    std::string sequence;          ///< Modified sequence (unmodified if the modified one is empty); empty for a non-peptide compound.
    std::vector<std::string> protein_refs;
  };

  /**
    @brief An OpenDIA spectral library: precursors and fragments as read from a PQP file.

    Only what OpenDIA uses is kept. Fragments are in the order of the file's transition query (explicitly sorted by
    precursor m/z, product m/z, RT, transition and precursor ID for libraries of at most 100,000 transitions).
  */
  struct Library
  {
    std::vector<LibraryCompound> compounds;
    std::vector<LibraryTransition> transitions;

    const std::vector<LibraryCompound>& getCompounds() const { return compounds; }
    const std::vector<LibraryTransition>& getTransitions() const { return transitions; }
  };

  /// Normalise a PQP ion type ("b-H20" is read as "b-H2O"; an unrecognised type becomes "unknown").
  std::string normaliseFragmentType(const std::string& type);

  /**
    @brief Read a PQP library (read-only SQLite access, the same query and row order as the OpenMS PQP reader).
    @throws Exception::FileNotReadable If the file cannot be opened.
    @throws Exception::InvalidValue If required tables are missing.
  */
  Library readLibrary(const std::string& pqp);

  /// PRECURSOR.ID to PRECURSOR.TRAML_ID of a PQP file.
  std::unordered_map<std::string, std::string> readPrecursorTraMLIDs(const std::string& pqp);
}
