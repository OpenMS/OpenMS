// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

namespace OpenMS
{
  /**
    @brief Tag names that TOPP tools attach to their parameters.

    The tags are stored in Param entries and interpreted by the parameter serializers
    (ParamCTDFile, ParamCWLFile) and by the tool framework. They live in the core
    library so that the serializers do not depend on TOPPBase, which is part of the
    separate OpenMS_CLI library; TOPPBase exposes the same values under its historical
    names (TOPPBase::TAG_INPUT_FILE etc.).
  */
  namespace ParamTags
  {
    inline constexpr const char* TAG_OUTPUT_FILE = "output file";
    inline constexpr const char* TAG_INPUT_FILE = "input file";
    inline constexpr const char* TAG_OUTPUT_DIR = "output dir";
    inline constexpr const char* TAG_OUTPUT_PREFIX = "output prefix";
    inline constexpr const char* TAG_ADVANCED = "advanced";
    inline constexpr const char* TAG_REQUIRED = "required";
  } // namespace ParamTags
} // namespace OpenMS
