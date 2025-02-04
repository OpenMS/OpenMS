// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
