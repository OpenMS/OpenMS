// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <set>

namespace OpenMS
{
    class FeatureMap;

    /**
      @brief Helper functions for FAIMS data

      FAIMSHelper contains convenience functions to deal with FAIMS
      compensation voltages and related data.

    */
    class OPENMS_DLLAPI FAIMSHelper
    {
    public:
      virtual ~FAIMSHelper() {}

      /**
        @brief Get all unique FAIMS compensation voltages (CVs) that occur in a PeakMap

        - Scans all spectra in the experiment and collects CVs from spectra whose drift time
          unit is DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE.
        - If the data does not contain any FAIMS spectra, an empty set will be returned.
        - The sentinel IMTypes::DRIFTTIME_NOT_SET is ignored; a warning is logged if encountered.

        @param[in] exp The input experiment
        @return Unique FAIMS compensation voltages (in volts)
      */
      static std::set<double> getCompensationVoltages(const PeakMap& exp);

      /**
        @brief Filter peptide identifications by FAIMS compensation voltage

        Filters peptide identifications to only include those matching the specified
        FAIMS CV. IDs without FAIMS_CV annotation are included for backward compatibility.

        @param[in] peptides Input peptide identifications
        @param[in] target_cv Target FAIMS compensation voltage to filter for
        @param[in] cv_tolerance Tolerance for floating point comparison (default: 0.01)
        @return Filtered list of peptide identifications
      */
      static PeptideIdentificationList filterPeptidesByFAIMSCV(
        const PeptideIdentificationList& peptides,
        double target_cv,
        double cv_tolerance = 0.01);

      /**
        @brief Filter identifications (identification data) by FAIMS compensation voltage

        As filterPeptidesByFAIMSCV(): keeps the identifications with the specified FAIMS CV (meta value FAIMS_CV) and
        those without it; the runs stay.

        @param[in] data Input identification data
        @param[in] target_cv Target FAIMS compensation voltage to filter for
        @param[in] cv_tolerance Tolerance for floating point comparison (default: 0.01)
        @return The filtered identification data
      */
      static IdentificationData filterIdentificationsByFAIMSCV(
        const IdentificationData& data,
        double target_cv,
        double cv_tolerance = 0.01);

      /**
        @brief Filter features by FAIMS compensation voltage

        Same rule as filterPeptidesByFAIMSCV: features whose FAIMS_CV matches @p target_cv are kept, and features
        without FAIMS_CV annotation are kept in every group. The document identifier and the identification data of
        @p features are copied (the kept features link it).

        @param[in] features Input features (e.g. untargeted seeds)
        @param[in] target_cv Target FAIMS compensation voltage to filter for
        @param[in] cv_tolerance Tolerance for floating point comparison (default: 0.01)
        @return Filtered features
      */
      static FeatureMap filterFeaturesByFAIMSCV(
        const FeatureMap& features,
        double target_cv,
        double cv_tolerance = 0.01);
    };

} //end namespace OpenMS
