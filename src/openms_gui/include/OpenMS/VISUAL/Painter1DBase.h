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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/VISUAL/PainterBase.h>


class QPainter;
class QPenStyle;

namespace OpenMS
{
  class LayerData1DBase;
  class LayerData1DChrom;
  class LayerData1DIonMobility;
  class LayerData1DPeak;
  class Plot1DCanvas;

  /**
   * @brief A base class for painting all items from a data layer (as supported by class derived from here) onto a 1D Canvas
  */
  class OPENMS_GUI_DLLAPI Painter1DBase : public PainterBase
  {
  public:
    virtual ~Painter1DBase() = default;

    /**
       @brief Paints items using the given painter onto the canvas.
 
       @param painter The painter used for drawing 
       @param canvas The canvas to paint onto (should expose all the details needed, like canvas size, draw mode, colors etc)
       @param layer_index Which layer is currently painted (FIXME: remove when Canvas1D::DrawMode and PenStyle are factored out) 
    */
    virtual void paint(QPainter* painter, Plot1DCanvas* canvas, int layer_index) = 0;

    void drawAnnotations_(const LayerData1DBase* layer, QPainter& painter, Plot1DCanvas* canvas) const;
  };

  /**
     @brief Painter1D for spectra
     
  */
  class OPENMS_GUI_DLLAPI Painter1DPeak : public Painter1DBase
  {
  public:
    /// C'tor which remembers the layer to paint
    Painter1DPeak(const LayerData1DPeak* parent);

    /// Implementation of base class
    void paint(QPainter*, Plot1DCanvas* canvas, int layer_index) override;

  protected:
    /// annotate up to 10 interesting peaks in the range @p vbegin to @p vend with their m/z values (using deisotoping and intensity filtering)
    void drawMZAtInterestingPeaks_(QPainter& painter, Plot1DCanvas* canvas, MSSpectrum::ConstIterator v_begin, MSSpectrum::ConstIterator v_end) const;

    const LayerData1DPeak* layer_; ///< the data to paint
  };

  /**
     @brief Painter1D for chromatograms

  */
  class OPENMS_GUI_DLLAPI Painter1DChrom : public Painter1DBase
  {
  public:
    /// C'tor which remembers the layer to paint
    Painter1DChrom(const LayerData1DChrom* parent);

    /// Implementation of base class
    void paint(QPainter*, Plot1DCanvas* canvas, int layer_index) override;

  protected:
    const LayerData1DChrom* layer_; ///< the data to paint
  };

    /**
   @brief Painter1D for mobilograms

*/
  class OPENMS_GUI_DLLAPI Painter1DIonMobility : public Painter1DBase
  {
  public:
    /// C'tor which remembers the layer to paint
    Painter1DIonMobility(const LayerData1DIonMobility* parent);

    /// Implementation of base class
    void paint(QPainter*, Plot1DCanvas* canvas, int layer_index) override;

  protected:
    const LayerData1DIonMobility* layer_; ///< the data to paint
  };

} // namespace OpenMS
