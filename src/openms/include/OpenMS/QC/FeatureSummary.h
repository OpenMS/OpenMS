// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/QC/QCBase.h>
#include <OpenMS/KERNEL/FeatureMap.h>

/**
 * @brief Detected Compounds as a Metabolomics QC metric
 *
 * Simple class to return a summary of detected compounds
 * from a featureXML file.
 *
 */

namespace OpenMS
{
  class OPENMS_DLLAPI FeatureSummary : public QCBase
  {
  public:
    /// Constructor
    FeatureSummary() = default;

    /// Destructor
    virtual ~FeatureSummary() = default;

    // stores feature summary values calculated by compute function
    struct OPENMS_DLLAPI Result
    {
      UInt feature_count = 0;
      float rt_shift_mean = 0;

      bool operator==(const Result& rhs) const;
    };

     /**
    @brief computes a summary of a featureXML file

    @param feature_map FeatureMap
    @return result object with summary values: 
            number of detected compounds (detected_compounds),
            retention time shift mean (rt_shift_mean)
    **/
    Result compute(const FeatureMap& feature_map);

    const String& getName() const override;

    QCBase::Status requires() const override;

  private:
    const String name_ = "Summary of features from featureXML file";
  };
}
