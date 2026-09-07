// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm, Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/FileTypes.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_MapAlignerBase MapAlignerBase

    @brief Base class for different MapAligner TOPP tools.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{
  // @brief stores model defaults for map aligner algorithms
  struct OPENMS_DLLAPI MapAlignerBase
  {
    static Param getModelDefaults(const std::string& default_model);
  };


class OPENMS_DLLAPI TOPPMapAlignerBase :
  public TOPPBase, public MapAlignerBase
{

public:
  TOPPMapAlignerBase(std::string name, std::string description, bool official = true) :
    TOPPBase(name, description, official), ref_params_(REF_NONE)
  {
  }

protected:

  // Kind of reference parameters that the tool offers:
  // - REF_NONE: no reference
  // - REF_RESTRICTED: reference file must have same type as input files
  // - REF_FLEXIBLE: reference file can have any supported file type
  enum ReferenceParameterKind { REF_NONE, REF_RESTRICTED, REF_FLEXIBLE }
    ref_params_;

  void registerOptionsAndFlagsMapAligners_(const std::string& file_formats,
                                           enum ReferenceParameterKind ref_params);

  ExitCodes checkParameters_();

  void transformSpectraFiles_(const StringList& in_spectra_files, 
                              const StringList& out_spectra_files,
                              const std::vector<TransformationDescription>& transformations,
                              bool store_original_rt);

};

}

/// @endcond

