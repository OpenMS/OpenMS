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

#include <OpenMS/VISUAL/LayerDataBase.h>

class QPainter;
class QPenStyle;

namespace OpenMS
{
  class LayerDataPeak;
  class Plot1DCanvas;

  /**
   * @brief A base class for painting all items from a data layer (as supported by class derived from here) onto a 1D Canvas
  */
  class OPENMS_GUI_DLLAPI Painter1DBase
  {
  public:
    virtual ~Painter1DBase() = default;

    /**
       @brief Paints items using the given painter onto the canvas.
 
       @param painter The painter used for drawing 
       @param canvas The canvas to paint onto (should expose all the details needed, like canvas size, draw mode, colors etc)
       @param layer_index Which layer is currently painted (FIXME: remove when Canvas1D::DrawMode and PenStyle are factored out) 
    */
    virtual void paint(QPainter*, Plot1DCanvas* canvas, int layer_index) = 0;

    /// static method to draw a dashed line
    static void drawDashedLine(const QPoint& from, const QPoint& to, QPainter* painter, QColor color);
  };

  /**
     @brief Painter1D for spectra
     
  */
  class OPENMS_GUI_DLLAPI Painter1DPeak : public Painter1DBase
  {
  public:
    /// C'tor which remembers the layer to paint
    Painter1DPeak(const LayerDataPeak* parent);

    /// Implementation of base class
    void paint(QPainter*, Plot1DCanvas* canvas, int layer_index) override;

  protected:
    /// draw all Annotation1DItems attached to the layer
    void drawAnnotations_(QPainter& painter, Plot1DCanvas* canvas);
    /// annotate up to 10 interesting peaks in the range @p vbegin to @pvend with their m/z values (using deisotoping and intensity filtering)
    void drawMZAtInterestingPeaks_(QPainter& painter, Plot1DCanvas* canvas, MSSpectrum::ConstIterator vbegin, MSSpectrum::ConstIterator vend);

    const LayerDataPeak* layer_; ///< the data to paint
  };

} // namespace OpenMS
