// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

/**
  * @brief FLASHIda C++ to C# (or vice versa) bridge functions
  * The functions here are called in C# to invoke functions in FLASHIda class
  * @see FLASHIda
  */

namespace OpenMS
{
  /// create FLASHIda class in C# FLASHIda side. Invoke FLASHIda constructor
  extern "C" OPENMS_DLLAPI FLASHIda *CreateFLASHIda(char *arg);

  /// delete FLASHIda class in C# FLASHIda side. Invoke FLASHIda destructor
  extern "C" OPENMS_DLLAPI void DisposeFLASHIda(FLASHIda *object);

  /// bridges getPeakGroups in FLASHIda class to C# FLASHIda side
  extern "C" OPENMS_DLLAPI int GetPeakGroupSize(FLASHIda *object,
                                                double *mzs,
                                                double *ints,
                                                int length,
                                                double rt_min,
                                                int ms_level,
                                                char *name,
                                                char *cv);

  /// bridges getIsolationWindows in FLASHIda class to C# FLASHIda side
  extern "C" OPENMS_DLLAPI void GetIsolationWindows(FLASHIda *object,
                                                    double *wstart,
                                                    double *wend,
                                                    double *qscores,
                                                    int *charges,
                                                    int *min_charges,
                                                    int *max_charges,
                                                    double *mono_masses,
                                                    double *chare_cos,
                                                    double *charge_snrs,
                                                    double *iso_cos,
                                                    double *snrs, double *charge_scores,
                                                    double *ppm_errors,
                                                    double *precursor_intensities,
                                                    double *peakgroup_intensities,
                                                    int *ids);

  // bridges removeFromExclusionList in FLASHida class to C# FLASHIda side
  extern "C" OPENMS_DLLAPI void RemoveFromExclusionList(FLASHIda *object, int id);

  extern "C" OPENMS_DLLAPI int GetAllPeakGroupSize(FLASHIda *pObject);

  extern "C" OPENMS_DLLAPI void GetAllMonoisotopicMasses(FLASHIda *pObject, double *masses, int length);

  extern "C" OPENMS_DLLAPI double GetRepresentativeMass(FLASHIda *pObject);

  /// keeps the precalculated averagine to calculate average masses from monoisotopic masses
  static FLASHHelperClasses::PrecalculatedAveragine avg;
}
