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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MRM_MRMFRAGMENTSELECTION_H
#define OPENMS_ANALYSIS_MRM_MRMFRAGMENTSELECTION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

namespace OpenMS
{
  /**
      @brief This class can select appropriate fragment ions of an MS/MS spectrum of a peptide

      @htmlinclude OpenMS_MRMFragmentSelection.parameters

      Several user choices can influence the selection of the ions from the MS/MS spectrum. These
      choices can be done using the parameters as described on the parameters page (see below).
      Basically there are two different ways of selecting suitable ions. One, using standardized
      names, e.g. given in the meta value "IonName" of each peaks of the spectrum (this can be
      written from TheoreticalSpectrumGenerator, PILISModel...). The second one is simply using
      the most abundant peaks in a specified m/z range.

      @ingroup Analysis_MRM
  */
  class OPENMS_DLLAPI MRMFragmentSelection :
    public DefaultParamHandler
  {

public:

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    MRMFragmentSelection();

    /// copy constructor
    MRMFragmentSelection(const MRMFragmentSelection & rhs);

    /// destructor
    ~MRMFragmentSelection() override;
    //@}

    /// assignment operator
    MRMFragmentSelection & operator=(const MRMFragmentSelection & rhs);

    /// selects accordingly to the parameters the best peaks of spec and writes them into selected_peaks
    void selectFragments(std::vector<Peak1D> & selected_peaks, const PeakSpectrum & spec);

protected:

    /// returns true if the selection of peak is allowed, according to the parameters set and the ion name
    bool peakselectionIsAllowed_(const String& name, const int charge);
  };
}

#endif
