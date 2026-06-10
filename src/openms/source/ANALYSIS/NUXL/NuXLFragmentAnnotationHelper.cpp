// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentAnnotationHelper.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <set>
#include <map>
#include <vector>
#include <algorithm>

using namespace OpenMS;


namespace OpenMS
{
  std::string NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon(char c, const std::string& fragment_shift_name)
  {
    return std::string("i") + c + "+" + fragment_shift_name + "+";
  }

  std::vector<PeptideHit::PeakAnnotation> NuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToPHFA(
    const std::string& ion_type, 
    const std::map<Size, std::vector<FragmentAnnotationDetail_> >& ion_annotation_details)
  {
    std::vector<PeptideHit::PeakAnnotation> fas;
    for (const auto& ait : ion_annotation_details)
    {
      for (const auto& sit : ait.second)
      {
        PeptideHit::PeakAnnotation fa;
        fa.charge = sit.charge;
        fa.mz = sit.mz;
        fa.intensity = sit.intensity;
        if (sit.shift.empty())
        {
          fa.annotation = ion_type + StringUtils::toStr(ait.first) + std::string(fa.charge, '+');
        }
        else
        {
          const std::string annotation_text = ion_type + StringUtils::toStr(ait.first) + "+" + sit.shift + std::string(fa.charge, '+'); 
          fa.annotation = annotation_text;
        }
        fas.push_back(std::move(fa));
      }
    }
    return fas;
  }

   std::vector<PeptideHit::PeakAnnotation> NuXLFragmentAnnotationHelper::shiftedToPHFA(
    const std::map<std::string, 
    std::set<std::pair<std::string, double> > >& shifted_ions)
  {
    std::vector<PeptideHit::PeakAnnotation> fas;
    for (const auto& ait : shifted_ions)
    {
      for (const auto& sit : ait.second)
      {
        PeptideHit::PeakAnnotation fa;
        fa.charge = 1;
        fa.mz = sit.second;
        fa.intensity = 1;
        const std::string annotation_text = sit.first;
        fa.annotation = annotation_text;
        fas.push_back(std::move(fa)); 
      }
    }
    return fas;
  }

  std::string NuXLFragmentAnnotationHelper::shiftedIonsToString(const std::vector<PeptideHit::PeakAnnotation>& as)
  {
    std::vector<PeptideHit::PeakAnnotation> sorted(as);
    stable_sort(sorted.begin(), sorted.end());
    std::string fas;
    for (const auto & a : sorted)
    {
      fas +=std::string("(") + StringUtils::number(a.mz, 3) + "," 
        + StringUtils::number(100.0 * a.intensity, 1) + ",\"" 
        + a.annotation + "\")";    
      if (&a != &sorted.back()) { fas += "|"; }     
    }
    return fas;
  }

  void NuXLFragmentAnnotationHelper::addShiftedPeakFragmentAnnotation_(
                                        const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_b_ions,
                                        const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_y_ions,
                                        const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_a_ions,
                                        const std::vector<PeptideHit::PeakAnnotation>& shifted_immonium_ions,
                                        const std::vector<PeptideHit::PeakAnnotation>& annotated_marker_ions,
                                        const std::vector<PeptideHit::PeakAnnotation>& annotated_precursor_ions,
                                        std::vector<PeptideHit::PeakAnnotation>& fas) 
  {
    if (!shifted_b_ions.empty())
    {
      const auto& fas_tmp = fragmentAnnotationDetailsToPHFA("b", shifted_b_ions);
      fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
    }

    if (!shifted_y_ions.empty())
    {
      const auto& fas_tmp = fragmentAnnotationDetailsToPHFA("y", shifted_y_ions);
      fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
    }

    if (!shifted_a_ions.empty())
    {
      const auto& fas_tmp = fragmentAnnotationDetailsToPHFA("a", shifted_a_ions);
      fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
    }

    if (!shifted_immonium_ions.empty())
    {
      fas.insert(fas.end(), shifted_immonium_ions.begin(), shifted_immonium_ions.end());
    }

    if (!annotated_marker_ions.empty())
    {
      fas.insert(fas.end(), annotated_marker_ions.begin(), annotated_marker_ions.end());
    }

    if (!annotated_precursor_ions.empty())
    {
      fas.insert(fas.end(), annotated_precursor_ions.begin(), annotated_precursor_ions.end());
    }
  }

} // namespace

