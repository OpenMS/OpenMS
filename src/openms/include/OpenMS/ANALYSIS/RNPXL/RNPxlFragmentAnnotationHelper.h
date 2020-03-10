// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
namespace Internal
{

/** 
   @brief Convenience functions to construct appealing fragment annotation strings
         and store them as PeptideHit::PeakAnnotation
 
 */
class OPENMS_DLLAPI RNPxlFragmentAnnotationHelper
{
  public:

  /// Single fragment annotation
  struct OPENMS_DLLAPI FragmentAnnotationDetail_
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

  static String getAnnotatedImmoniumIon(char c, const String& fragment_shift_name);

  /// conversion of RNPxl annotations to PeptideHit::PeakAnnotation
  static std::vector<PeptideHit::PeakAnnotation> fragmentAnnotationDetailsToPHFA(
    const String& ion_type, 
    std::map<Size, std::vector<FragmentAnnotationDetail_> > ion_annotation_details);

  static std::vector<PeptideHit::PeakAnnotation> shiftedToPHFA(
    const std::map<String, 
    std::set<std::pair<String, double> > >& shifted_ions);


  static String shiftedIonsToString(const std::vector<PeptideHit::PeakAnnotation>& as);

  static void addShiftedPeakFragmentAnnotation_(const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_b_ions,
                                         const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_y_ions,
                                         const std::map<Size, std::vector<FragmentAnnotationDetail_>>& shifted_a_ions,
                                         const std::vector<PeptideHit::PeakAnnotation>& shifted_immonium_ions,
                                         const std::vector<PeptideHit::PeakAnnotation>& annotated_marker_ions,
                                         const std::vector<PeptideHit::PeakAnnotation>& annotated_precursor_ions,
                                         std::vector<PeptideHit::PeakAnnotation>& fas);
};
} // namespace Internal
} // namespace OpenMS


