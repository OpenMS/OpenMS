// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <filesystem>
#include <string>
#include <vector>

/**
  Private helpers shared by FeatureMapArrowIO and ConsensusMapArrowIO.

  A feature or consensus bundle stores the map's owning identification data as a nested
  native bundle (identifications/, see IdentificationDataFile) and the per-feature links
  (primary molecule, observed queries, matches) in identification_links.parquet, keyed by
  feature unique ID like psms.parquet. Both are optional: bundles without native
  identifications, including those written by OpenMS 3.6, contain neither.
*/
namespace OpenMS::Internal::MapIdentificationParquet
{
  inline const std::string IDENTIFICATIONS_DIRECTORY = "identifications";
  inline const std::string LINKS_FILE = "identification_links.parquet";

  /// True if the map owns identification data or any feature carries a native link.
  bool hasNativeIdentifications(const IdentificationData& data, const std::vector<const BaseFeature*>& features);

  /// Throws Exception::InvalidValue unless every query and match link resolves in @p data.
  void validateLinks(const IdentificationData& data, const std::vector<const BaseFeature*>& features);

  /**
    Writes identifications/ and identification_links.parquet into @p directory (a fresh staging
    directory). Linked features need distinct, valid unique IDs. Throws on failure.
  */
  void store(const std::filesystem::path& directory, const IdentificationData& data, const std::vector<const BaseFeature*>& features);

  /**
    Loads identifications/ and identification_links.parquet from @p directory if present and
    applies the links to @p features by unique ID. Leaves @p data and @p features unchanged on failure.
  */
  void load(const std::filesystem::path& directory, IdentificationData& data, const std::vector<BaseFeature*>& features);

  /**
    Writes a map bundle through a temporary sibling directory and publishes it atomically.

    The target may be absent, an empty directory or an existing bundle of the same kind (identified
    by @p main_file), which is replaced; any other existing file or directory is rejected.
  */
  class StagedBundle
  {
  public:
    StagedBundle(const std::string& target, const std::string& main_file);
    ~StagedBundle();
    StagedBundle(const StagedBundle&) = delete;
    StagedBundle& operator=(const StagedBundle&) = delete;
    const std::filesystem::path& path() const
    { return path_; }
    void publish();

  private:
    std::filesystem::path sibling_(const std::string& infix) const;
    bool replaceable_() const;
    std::filesystem::path target_;
    std::string main_file_;
    std::filesystem::path path_;
    bool published_ = false;
  };
} // namespace OpenMS::Internal::MapIdentificationParquet
