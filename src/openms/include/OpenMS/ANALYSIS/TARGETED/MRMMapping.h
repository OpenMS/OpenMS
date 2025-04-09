// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  /**
    @brief A class to map targeted assays to chromatograms.

  */
  class OPENMS_DLLAPI MRMMapping :
    public DefaultParamHandler
  {
public:

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    MRMMapping();

    /// destructor
    ~MRMMapping() override {}
    //@}

    /**
      @brief Maps input chromatograms to assays in a targeted experiment
      
      The output chromatograms are an annotated copy of the input chromatograms
      with native id, precursor information and peptide sequence (if available)
      annotated in the chromatogram files.

      The algorithm tries to match a given set of chromatograms and targeted
      assays. It iterates through all the chromatograms retrieves one or more
      matching targeted assay for the chromatogram. By default, the algorithm
      assumes that a 1:1 mapping exists. If a chromatogram cannot be mapped
      (does not have a corresponding assay) the algorithm issues a warning, the
      user can specify that the program should abort in such a case (see
      error_on_unmapped).
      
      @note If multiple mapping is enabled (see map_multiple_assays parameter)
      then each mapped assay will get its own chromatogram that contains the
      same raw data but different meta-annotation. This *can* be useful if the
      same transition is used to monitor multiple analytes but may also
      indicate a problem with too wide mapping tolerances.
    */
    void mapExperiment(const OpenMS::PeakMap& input_chromatograms,
        const OpenMS::TargetedExperiment& targeted_exp,
        OpenMS::PeakMap& output) const;

protected:

    /// copy constructor
    MRMMapping(const MRMMapping & rhs);

    /// assignment operator
    MRMMapping & operator=(const MRMMapping & rhs);

    /// Synchronize members with param class
    void updateMembers_() override;

    double precursor_tol_;
    double product_tol_;
    bool map_multiple_assays_;
    bool error_on_unmapped_;

  };
}

