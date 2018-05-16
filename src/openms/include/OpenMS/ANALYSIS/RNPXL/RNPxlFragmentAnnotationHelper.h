// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_RNPXL_FRAGMENTANNOTATIONHELPER_H
#define OPENMS_ANALYSIS_RNPXL_FRAGMENTANNOTATIONHELPER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <set>
#include <map>
#include <vector>
#include <algorithm>

namespace OpenMS 
{

/* @brief Convenience functions to construct appealing fragment annotation strings
 *
 */
class FragmentAnnotationHelper
{
  public:

  /// Single fragment annotation
  struct FragmentAnnotationDetail_
  {
    FragmentAnnotationDetail_(String s, int z, double m, double i):
      shift(s),
      charge(z),
      mz(m),
      intensity(i)
      {}
    String shift;
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

  static String getAnnotatedImmoniumIon(char c, const String& fragment_shift_name)
  {
    return String("i") + c + "+" + fragment_shift_name;
  }

  // deprecated: for PD community nodes compatibility 
  static String fragmentAnnotationDetailsToString(const String& ion_type, std::map<Size, std::vector<FragmentAnnotationDetail_> > ion_annotation_details)
  {
    String fas;
    for (auto ait = ion_annotation_details.begin(); ait != ion_annotation_details.end(); ++ait)
    {
      for (auto sit = ait->second.begin(); sit != ait->second.end(); ++sit)
      {
        if (ait != ion_annotation_details.begin() || sit != ait->second.begin())
        {
          fas += "|";
        }

        String annotation_text;
        annotation_text = sit->shift.empty() ? "[" + ion_type + String(ait->first) + "]" + String(sit->charge, '+') : "[" + ion_type + String(ait->first) + "+" + sit->shift + "]" + String(sit->charge, '+'); // e.g.: [b3]+ and  [y3+H3PO4]++
        // e.g.: (343.5,99.5,"[b2-H2O]+")
        fas += "(" + String::number(sit->mz, 3) + "," + String::number(100.0 * sit->intensity, 1) + "," + "\"" + annotation_text+ "\")";
      }
    }
    return fas;
  }

  // conversion of RNPxl annotations to PeptideHit::PeakAnnotation
  static std::vector<PeptideHit::PeakAnnotation> fragmentAnnotationDetailsToPHFA(const String& ion_type, std::map<Size, std::vector<FragmentAnnotationDetail_> > ion_annotation_details)
  {
    std::vector<PeptideHit::PeakAnnotation> fas;
    for (auto ait : ion_annotation_details)
    {
      for (auto sit : ait.second)
      {
        PeptideHit::PeakAnnotation fa;
        fa.charge = sit.charge;
        fa.mz = sit.mz;
        fa.intensity = sit.intensity;
        if (sit.shift.empty())
        {
          fa.annotation = ion_type + String(ait.first);
        }
        else
        {
          const String annotation_text = ion_type + String(ait.first) + "+" + sit.shift; 
          fa.annotation = annotation_text;
        }
        fas.push_back(fa);
      }
    }
    return fas;
  }

  static std::vector<PeptideHit::PeakAnnotation> shiftedToPHFA(const std::map<String, std::set<std::pair<String, double> > >& shifted_ions)
  {
    std::vector<PeptideHit::PeakAnnotation> fas;
    for (auto ait : shifted_ions)
    {
      for (auto sit : ait.second)
      {
        PeptideHit::PeakAnnotation fa;
        fa.charge = 1;
        fa.mz = sit.second;
        fa.intensity = 1;
        const String annotation_text = sit.first;
        fa.annotation = annotation_text;
        fas.push_back(fa); 
      }
    }
    return fas;
  }


  static String shiftedIonsToString(const std::vector<PeptideHit::PeakAnnotation>& as)
  {
    std::vector<PeptideHit::PeakAnnotation> sorted(as);
    stable_sort(sorted.begin(), sorted.end());
    String fas;
    for (auto & a : sorted)
    {
      fas += String("(") + String::number(a.mz, 3) + "," + String::number(100.0 * a.intensity, 1) + ",\"" + a.annotation + "\")";    
      if (&a != &sorted.back()) { fas += "|"; }     
    }
    return fas;
  }

  static void addShiftedPeakFragmentAnnotation_(const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_b_ions,
                                         const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_y_ions,
                                         const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_a_ions,
                                         const std::vector<PeptideHit::PeakAnnotation>& shifted_immonium_ions,
                                         const std::vector<PeptideHit::PeakAnnotation>& annotated_marker_ions,
                                         const std::vector<PeptideHit::PeakAnnotation>& annotated_precursor_ions,
                                         StringList &fa_strings, 
                                         std::vector<PeptideHit::PeakAnnotation>& fas) 
  {
    String sb = fragmentAnnotationDetailsToString("b", shifted_b_ions);
    if (!sb.empty())
    {
      fa_strings.push_back(sb);
      const std::vector<PeptideHit::PeakAnnotation>& fas_tmp = fragmentAnnotationDetailsToPHFA("b", shifted_b_ions);;
      fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
    }

    String sy = fragmentAnnotationDetailsToString("y", shifted_y_ions);
    if (!sy.empty())
    {
      fa_strings.push_back(sy);
      const std::vector<PeptideHit::PeakAnnotation>& fas_tmp = fragmentAnnotationDetailsToPHFA("y", shifted_y_ions);;
      fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
    }

    String sa = fragmentAnnotationDetailsToString("a", shifted_a_ions);
    if (!sa.empty())
    {
      fa_strings.push_back(sa);
      const std::vector<PeptideHit::PeakAnnotation>& fas_tmp = fragmentAnnotationDetailsToPHFA("a", shifted_a_ions);;
      fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
    }

    String sii = shiftedIonsToString(shifted_immonium_ions);
    if (!sii.empty())
    {
      fa_strings.push_back(sii);
      fas.insert(fas.end(), shifted_immonium_ions.begin(), shifted_immonium_ions.end());
    }

    String smi = shiftedIonsToString(annotated_marker_ions);
    if (!smi.empty())
    {
      fa_strings.push_back(smi);
      fas.insert(fas.end(), annotated_marker_ions.begin(), annotated_marker_ions.end());
    }

    String spi = shiftedIonsToString(annotated_precursor_ions);
    if (!spi.empty())
    {
      fa_strings.push_back(spi);
      fas.insert(fas.end(), annotated_precursor_ions.begin(), annotated_precursor_ions.end());
    }
  }
};

}

#endif

