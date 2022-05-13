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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include "Plot1DCanvas.h"

#include <OpenMS/VISUAL/LayerDataBase.h>

class QWidget;

namespace OpenMS
{
  class Annotations1DContainer;
  class Annotation1DItem;

  /**
   * \brief Base class for all 1D layers, a special case of LayerData.
   *
   * 1D is a bit special because we need to remember which spectrum/chrom/IM is currently shown (there are usually many of them to choose from).
   *
   *
   */
  class OPENMS_DLLAPI LayerData1DBase : public virtual LayerDataBase
  {
  public:
    LayerData1DBase() {};

    virtual std::unique_ptr<Painter1DBase> getPainter1D() const = 0;

    /// get name augmented with attributes, e.g. '*' if modified
    String getDecoratedName() const override;
        
    /// Returns a const reference to the annotations of the current spectrum (1D view)
    const Annotations1DContainer& getCurrentAnnotations() const
    {
      return annotations_1d_[current_idx_];
    }

    /// Returns a mutable reference to the annotations of the current spectrum (1D view)
    Annotations1DContainer& getCurrentAnnotations()
    {
      return annotations_1d_[current_idx_];
    }

    /// Returns a const reference to the annotations of the @p spectrum_index's spectrum (1D view)
    const Annotations1DContainer& getAnnotations(Size spectrum_index) const
    {
      return annotations_1d_[spectrum_index];
    }

    /// Returns a mutable reference to the annotations of the @p spectrum_index's spectrum (1D view)
    Annotations1DContainer& getAnnotations(Size spectrum_index)
    {
      return annotations_1d_[spectrum_index];
    }
    /// Get the index of the current spectrum (1D view)
    Size getCurrentIndex() const
    {
      return current_idx_;
    }

    /// Set the index of the current spectrum (1D view)
    void setCurrentIndex(Size index)
    {
      current_idx_ = index;
    }


    /// if this layer is flipped (1d mirror view)
    bool flipped = false;

    /// Peak colors of the currently shown spectrum
    std::vector<QColor> peak_colors_1d;

  protected:
    /// Index of the current spectrum/chromatogram etc
    Size current_idx_ = -1;

    /// Annotations of all spectra of the experiment (1D view)
    std::vector<Annotations1DContainer> annotations_1d_ = std::vector<Annotations1DContainer>(1);


  };


}// namespace OpenMS

