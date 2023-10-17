// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

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

