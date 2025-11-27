// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerData1DBase.h>

#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
                                                        
using namespace std;

namespace OpenMS
{
  String LayerData1DBase::getDecoratedName() const
  {
    String n = LayerDataBase::getDecoratedName();
    if (flipped)
    {
      n += " [flipped]";
    }
    return n;
  }

  void LayerData1DBase::setCurrentIndex(Size index)
  {
    current_idx_ = index;
    if (current_idx_ >= annotations_1d_.size())
    {
      annotations_1d_.resize(current_idx_ + 1);
    }
    
    // Clear peak colors to force reinitialization for the new spectrum
    // Unlike annotations which persist across spectra, peak colors need to be regenerated
    // to match the size of the new spectrum, preventing "Peak color array size doesn't
    // match number of peaks" errors that occur when switching between spectra with
    // different numbers of peaks
    peak_colors_1d.clear();
  }
}// namespace OpenMS
