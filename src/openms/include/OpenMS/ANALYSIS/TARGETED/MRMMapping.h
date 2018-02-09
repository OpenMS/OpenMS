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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_TARGETED_MRMMAPPING_H
#define OPENMS_ANALYSIS_TARGETED_MRMMAPPING_H

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
      
      @Note If multiple mapping is enabled (see map_multiple_assays parameter)
      then each mapped assay will get its own chromatogram that contains the
      same raw data but different meta-annotation. This *can* be useful if the
      same transition is used to monitor multiple analytes but may also
      indicate a problem with too wide mapping tolerances.
    */
    void mapExperiment(const OpenMS::PeakMap& input_chromatograms,
        const OpenMS::TargetedExperiment& targeted_exp,
        OpenMS::PeakMap& output);

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

#endif // OPENMS_ANALYSIS_TARGETED_MRMMAPPING_H
