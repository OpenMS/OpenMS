// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/OpenMSConfig.h>

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace OpenMS::Internal
{
/// One "trailer extra" label/value pair of a Thermo scan.
struct OPENMS_DLLAPI ThermoTrailerEntry
{
  std::string label; ///< e.g. "Monoisotopic M/Z:"
  std::string value; ///< as reported by the instrument, numbers included
};

/// One reaction (isolation plus activation step) of a Thermo MSn scan event.
struct OPENMS_DLLAPI ThermoReaction
{
  std::optional<double> precursor_mass;   ///< isolation target m/z
  std::optional<double> isolation_width;  ///< isolation window width (Th)
  std::optional<double> isolation_offset; ///< isolation window offset from the target (Th)
  std::string activation;                 ///< vendor activation type name, e.g. "HigherEnergyCollisionalDissociation"
  std::optional<double> collision_energy; ///< vendor value; meaningful only if collision_energy_valid
  bool collision_energy_valid = false;    ///< vendor flag for collision_energy
};

/// The per-scan metadata needed to reconstruct the precursor hierarchy.
struct OPENMS_DLLAPI ThermoScan
{
  int scan_number = 0;                     ///< 1-based scan number within the controller
  int ms_level = 0;                        ///< MS order (1 for full scans)
  std::string filter;                      ///< scan filter string, e.g. "FTMS + c NSI Full ms2 500.00@hcd30.00 [100-1000]"
  std::string native_id;                   ///< e.g. "controllerType=0 controllerNumber=1 scan=12"
  std::vector<ThermoTrailerEntry> trailer; ///< trailer-extra entries in file order
  std::vector<ThermoReaction> reactions;   ///< reactions in acquisition order (MS2 first)
};

/// One precursor of an MSn scan as reconstructed by ThermoRawFileMetadata::precursors().
struct OPENMS_DLLAPI ThermoPrecursor
{
  double target_mz = 0.0;         ///< isolation window target m/z
  double selected_mz = 0.0;       ///< monoisotopic m/z when plausible, otherwise the target
  std::optional<int> charge;      ///< unset if unknown or not positive
  std::optional<double> width;    ///< isolation width; unset if unknown or negative
  double lower_offset = 0.0;      ///< width / 2 - isolation offset (0 if the width is unknown)
  double upper_offset = 0.0;      ///< width / 2 + isolation offset
  int parent_scan = 0;            ///< scan this precursor was isolated from (0 if unknown)
  bool estimate_intensity = true; ///< false for additional SPS selections
  std::string spectrum_ref;       ///< native ID of the parent scan ("" if unknown)
  ThermoReaction activation;      ///< primary reaction
  std::optional<ThermoReaction> supplemental; ///< supplemental activation (e.g. the HCD step of EThcD)
};

/** @brief Resolve the precursor hierarchy from bridge metadata without fetching
   peaks. Keeps isolation targets distinct from selected ions, and retains SPS
   selections and supplemental reactions. The filter fallback follows
   ThermoRawFileParser. */
class OPENMS_DLLAPI ThermoRawFileMetadata
{
public:
  /// Value of the trailer-extra entry with exactly this @p label, or "" if absent.
  static std::string trailer(const ThermoScan& scan, const std::string& label);

  /**
    @brief Parse a trailer value as a number.
    @return the number, or an empty optional when the value is empty, not fully numeric (e.g. "2,5"
            or "12.3abc") or not finite. A literal zero is a valid number, not "missing".
  */
  static std::optional<double> number(const std::string& value);

  /**
    @brief Selected ion m/z of a precursor (ThermoRawFileParser rule).
    @param target isolation window target m/z of the reaction
    @param mono "Monoisotopic M/Z:" trailer value as returned by number()
    @param width isolation width (0 if unknown)
    @return the monoisotopic m/z when it is set and plausibly inside the isolation window
            (narrow windows accept -3.0 / +2.5 Th around the target), otherwise the target
  */
  static double selectedIon(double target, std::optional<double> mono, double width);

  /**
    @brief Precursor descriptors of one scan; scans must be passed in acquisition order.

    MS1 scans return an empty vector and are only recorded as potential parents. For MSn
    scans the vector holds the scan's own precursor (isolation target, selected ion, charge,
    window, parent scan / spectrum reference, activation and supplemental activation),
    followed by further SPS selections and then the precursors of every ancestor scan, so
    an MS3 scan lists its MS2 parent's precursor as well. Only each scan's own precursors
    are retained internally; ancestors are resolved through parent links on request.
  */
  std::vector<ThermoPrecursor> precursors(const ThermoScan& scan);

private:
  /// What is retained per seen scan to resolve later scans' parents.
  struct State
  {
    int level = 1;             ///< MS order
    int parent = 0;            ///< scan this one descends from (0 if unknown)
    std::size_t reactions = 0; ///< reactions consumed up to and including this scan
    std::string native_id;     ///< spectrum reference for descendants
    std::vector<ThermoPrecursor> own; ///< this scan's own precursor descriptors
  };
  std::map<int, State> scans_;         ///< seen scans by scan number
  std::map<std::string, int> filters_; ///< most recent scan per filter key ("" for MS1)

  /// True if @p second is the supplemental (HCD / CID) step of an ETD / ECD reaction @p first on the same target.
  static bool supplemental_(const ThermoReaction& first, const ThermoReaction& second);
  /// Index of the reaction that isolated this scan's precursor when no parent is known.
  static std::size_t lastReaction_(const std::vector<ThermoReaction>& reactions);
  /// Most recent scan whose filter key is the parent of @p key (0 if none).
  int parentFromFilter_(const std::string& key) const;
  /// SPS target masses from the legacy ("SPS Mass N:") or modern ("SPS Masses:") trailer labels.
  static std::vector<double> spsMasses_(const ThermoScan& scan);
};
} // namespace OpenMS::Internal
