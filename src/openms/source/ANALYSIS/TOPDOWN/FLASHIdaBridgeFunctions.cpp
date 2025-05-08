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
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIdaBridgeFunctions.h>

namespace OpenMS
{
  FLASHIda *CreateFLASHIda(char *arg)
  {
    std::cout << " FLASHIda creating ... " << std::endl;
    return new FLASHIda(arg);
  }

  void DisposeFLASHIda(FLASHIda *pObject)
  {
    if (pObject != nullptr)
    {
      delete pObject;
      pObject = nullptr;
    }
  }

  int GetPeakGroupSize(FLASHIda *pObject, double *mzs, double *ints, int length, double rt_min, int msLevel, char *name, char *cv)
  {
    if (pObject != nullptr)
    {
      return pObject->getPeakGroups(mzs, ints, length, rt_min * 60.0, msLevel, name, cv);
    }
    return 0;
  }

  double GetRepresentativeMass(FLASHIda *pObject)
  {
    return pObject->getRepresentativeMass();
  }

  void GetIsolationWindows(FLASHIda *pObject,
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
                           double *snrs,
                           double *charge_scores,
                           double *ppm_errors,
                           double *precursor_intensities,
                           double *peakgroup_intensities,
                           int *ids
  )
  {
    if (pObject != nullptr)
    {
      pObject->getIsolationWindows(wstart,
                                   wend,
                                   qscores,
                                   charges,
                                   min_charges,
                                   max_charges,
                                   mono_masses, chare_cos,
                                   charge_snrs, iso_cos, snrs, charge_scores, ppm_errors,
                                   precursor_intensities,
                                   peakgroup_intensities,
                                   ids);
    }
  }

  void RemoveFromExclusionList(
      FLASHIda *pObject,
      int id
  )
  {
      if (pObject != nullptr) 
      {
        pObject->removeFromExlusionList(id);
      }
  }

  int GetAllPeakGroupSize(FLASHIda *pObject)
  {if (pObject != nullptr)
    {
      return pObject->GetAllPeakGroupSize();
    }
    return 0;
  }

  void GetAllMonoisotopicMasses(FLASHIda *pObject, double *masses, int length)
  {
    if (pObject != nullptr)
    {
      pObject->getAllMonoisotopicMasses(masses, length);
    }
  }
}
