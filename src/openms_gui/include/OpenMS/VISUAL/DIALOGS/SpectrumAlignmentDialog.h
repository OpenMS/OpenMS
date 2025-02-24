// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <QtWidgets/QDialog>

namespace Ui
{
  class SpectrumAlignmentDialogTemplate;
}

namespace OpenMS
{
  class Plot1DWidget;

  /**
      @brief Lets the user select two spectra and set the parameters for the spectrum alignment.

      @ingroup Dialogs
  */
  class SpectrumAlignmentDialog :
    public QDialog
  {
    Q_OBJECT

public:

    /// Constructor
    SpectrumAlignmentDialog(Plot1DWidget * parent);
    ~SpectrumAlignmentDialog() override;
  
    double getTolerance() const;
    bool isPPM() const;

    /// Returns the index of the selected non-flipped layer
    Int get1stLayerIndex();
    /// Returns the index of the selected flipped layer
    Int get2ndLayerIndex();

protected slots:

protected:

    /// Stores the layer indices of the layers in the left list (non-flipped layers)
    std::vector<UInt> layer_indices_1_;
    /// Stores the layer indices of the layers in the right list (flipped layers)
    std::vector<UInt> layer_indices_2_;

private:
    Ui::SpectrumAlignmentDialogTemplate* ui_;
  };

}
