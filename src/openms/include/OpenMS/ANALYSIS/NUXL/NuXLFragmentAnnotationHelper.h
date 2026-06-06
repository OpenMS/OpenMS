// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <set>
#include <map>
#include <vector>
#include <algorithm>

namespace OpenMS
{

/** 
   @brief Convenience functions to construct appealing fragment annotation strings
         and store them as PeptideHit::PeakAnnotation
 
 */
class OPENMS_DLLAPI NuXLFragmentAnnotationHelper
{
  public:

  /// Single fragment annotation
  struct OPENMS_DLLAPI FragmentAnnotationDetail_
  {
    FragmentAnnotationDetail_(std::string s, int z, double m, double i):
      shift(s),
      charge(z),
      mz(m),
      intensity(i)
      {}
    std::string shift;
    int charge;
    double mz;
    double intensity;

    bool operator<(const FragmentAnnotationDetail_& other) const
    {
      return std::tie(charge, shift, mz, intensity) < 
        std::tie(other.charge, other.shift, other.mz, other.intensity);
    }

    bool operator==(const FragmentAnnotationDetail_& other) const
    {
      double mz_diff = fabs(mz - other.mz);
      double intensity_diff = fabs(intensity - other.intensity);
      return (charge == other.charge && shift == other.shift && mz_diff < 1e-6 && intensity_diff < 1e-6); // mz and intensity difference comparison actually not needed but kept for completeness
    }
  };

  static std::string getAnnotatedImmoniumIon(char c, const std::string& fragment_shift_name);

  /// conversion of NuXL annotations to PeptideHit::PeakAnnotation
  static std::vector<PeptideHit::PeakAnnotation> fragmentAnnotationDetailsToPHFA(
    const std::string& ion_type, 
    const std::map<Size, std::vector<FragmentAnnotationDetail_> >& ion_annotation_details);

  static std::vector<PeptideHit::PeakAnnotation> shiftedToPHFA(
    const std::map<std::string, 
    std::set<std::pair<std::string, double> > >& shifted_ions);


  static std::string shiftedIonsToString(const std::vector<PeptideHit::PeakAnnotation>& as);

  static void addShiftedPeakFragmentAnnotation_(const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_b_ions,
                                         const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_y_ions,
                                         const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_a_ions,
                                         const std::vector<PeptideHit::PeakAnnotation>& shifted_immonium_ions,
                                         const std::vector<PeptideHit::PeakAnnotation>& annotated_marker_ions,
                                         const std::vector<PeptideHit::PeakAnnotation>& annotated_precursor_ions,
                                         std::vector<PeptideHit::PeakAnnotation>& fas);
};
} // namespace OpenMS


