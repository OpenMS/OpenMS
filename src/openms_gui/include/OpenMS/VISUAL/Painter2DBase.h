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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/VISUAL/INTERFACES/IPeptideIds.h>
#include <OpenMS/VISUAL/PainterBase.h>

#include <vector>

class QPainter;
class QPenStyle;
class QPoint;
class String;

namespace OpenMS
{
  class ConsensusFeature;
  class LayerDataChrom;
  class LayerDataConsensus;
  class LayerDataFeature;
  class LayerDataIdent;
  class LayerDataIonMobility;
  class LayerDataPeak;
  struct PeakIndex;
  class Plot2DCanvas;


  /**
   * @brief A base class for painting all items from a data layer (as supported by class derived from here) onto a 2D Canvas
  */
  class OPENMS_GUI_DLLAPI Painter2DBase : public PainterBase
  {
  public:
    virtual ~Painter2DBase() = default;

    /**
       @brief Paints items using the given painter onto the canvas.
 
       @param painter The painter used for drawing 
       @param canvas The canvas to paint onto (should expose all the details needed, like canvas size, draw mode, colors etc)
       @param layer_index Which layer is currently painted
    */
    virtual void paint(QPainter* painter, Plot2DCanvas* canvas, int layer_index) = 0;

    /**
     * \brief Emphasize a certain element (e.g. feature), e.g. when mouse hovering.
     * By default, nothing is highlighted. Override for subclasses if you need highlighting.
     *
     * \param painter The painter used for drawing 
     * \param canvas The canvas to paint onto (should expose all the details needed, like canvas size, draw mode, colors etc)
     * \param element Which item of the current layer should be drawn?
     */
    virtual void highlightElement(QPainter* painter, Plot2DCanvas* canvas, const PeakIndex element);

  protected:
    /**
      @brief Paints a convex hull.

      @param painter The QPainter to paint on
      @param canvas The canvas (for configuration details)
      @param hull Reference to convex hull
      @param has_identifications Draw hulls in green (true) or blue color (false)
    */
    static void paintConvexHull_(QPainter& painter, Plot2DCanvas* canvas, const ConvexHull2D& hull, bool has_identifications);

    /**
      @brief Paints convex hulls.

      @param painter The QPainter to paint on
      @param canvas The canvas (for configuration details)
      @param hulls Reference to convex hulls
      @param has_identifications Draw hulls in green (true) or blue color (false)
    */
    static void paintConvexHulls_(QPainter& painter, Plot2DCanvas* canvas, const std::vector<ConvexHull2D>& hulls, bool has_identifications);

    static void paintPeptideIDs_(QPainter* painter, Plot2DCanvas* canvas, const IPeptideIds::PepIds& ids, int layer_index);
  };

  /**
     @brief Painter2D for spectra
     
  */
  class OPENMS_GUI_DLLAPI Painter2DPeak : public Painter2DBase
  {
  public:
    /// C'tor which remembers the layer to paint
    Painter2DPeak(const LayerDataPeak* parent);

    void paint(QPainter*, Plot2DCanvas* canvas, int layer_index) override;

  protected:
    void paintAllIntensities_(QPainter& painter, Plot2DCanvas* canvas, Size layer_index, double pen_width);

    /**
      @brief Paints maximum intensity of individual peaks.

      Paints the peaks as small ellipses. The peaks are colored according to the
      selected dot gradient.

      @param painter The QPainter to paint with.
      @param canvas The canvas to paint on.
      @param layer_index The index of the layer.
      @param rt_pixel_count
      @param mz_pixel_count
    */
    void paintMaximumIntensities_(QPainter& painter, Plot2DCanvas* canvas, Size layer_index, Size rt_pixel_count, Size mz_pixel_count);


    /**
      @brief Paints the locations where MS2 scans where triggered
    */
    void paintPrecursorPeaks_(QPainter& painter, Plot2DCanvas* canvas);
    const LayerDataPeak* layer_; ///< the data to paint
  };

  /**
     @brief Painter2D for chromatograms

  */
  class OPENMS_GUI_DLLAPI Painter2DChrom : public Painter2DBase
  {
  public:
    /// C'tor which remembers the layer to paint
    Painter2DChrom(const LayerDataChrom* parent);

    void paint(QPainter* painter, Plot2DCanvas* canvas, int layer_index) override;

  protected:
    const LayerDataChrom* layer_; ///< the data to paint
  };

 /**
   @brief Painter2D for ion mobilograms

  */
  class OPENMS_GUI_DLLAPI Painter2DIonMobility : public Painter2DBase
  {
  public:
    /// C'tor which remembers the layer to paint
    Painter2DIonMobility(const LayerDataIonMobility* parent);

    void paint(QPainter* painter, Plot2DCanvas* canvas, int layer_index) override;

  protected:
    const LayerDataIonMobility* layer_; ///< the data to paint
  };

  /**
     @brief Painter2D for Features

  */
  class OPENMS_GUI_DLLAPI Painter2DFeature : public Painter2DBase
  {
  public:
    /// C'tor which remembers the layer to paint
    Painter2DFeature(const LayerDataFeature* parent);

    void paint(QPainter*, Plot2DCanvas* canvas, int layer_index) override;

    void highlightElement(QPainter* painter, Plot2DCanvas* canvas, const PeakIndex element) override;

  protected:
    /**
      @brief Paints convex hulls (one for each mass trace) of a features layer.
    */
    void paintTraceConvexHulls_(QPainter* painter, Plot2DCanvas* canvas);

    /**
      @brief Paints the convex hulls (one for each feature) of a features layer.
    */
    void paintFeatureConvexHulls_(QPainter* painter, Plot2DCanvas* canvas);

    const LayerDataFeature* layer_; ///< the data to paint
  };

  /**
     @brief Painter2D for ConsensusFeatures

  */
  class OPENMS_GUI_DLLAPI Painter2DConsensus : public Painter2DBase
  {
  public:
    /// C'tor which remembers the layer to paint
    Painter2DConsensus(const LayerDataConsensus* parent);

    void paint(QPainter*, Plot2DCanvas* canvas, int layer_index) override;

    void highlightElement(QPainter* painter, Plot2DCanvas* canvas, const PeakIndex element) override;

  protected:
    /**
      @brief Paints the consensus elements of a consensus features layer.

      @param painter The QPainter to paint on
      @param canvas The canvas (for configuration details)
      @param layer_index Index of the layer
    */
    void paintConsensusElements_(QPainter* painter, Plot2DCanvas* canvas, Size layer_index);

    /**
      @brief Paints one consensus element of a consensus features layer.

      @param painter The QPainter to paint on
      @param canvas The canvas (for configuration details)
      @param layer_index Index of the layer
      @param cf Reference to the consensus feature to be painted
    */
    void paintConsensusElement_(QPainter* painter, Plot2DCanvas* canvas, Size layer_index, const ConsensusFeature& cf);

    /**
      @brief checks if any element of a consensus feature is currently visible.

      @param canvas The canvas (for configuration details)
      @param cf The ConsensusFeature that needs checking
      @param layer_index Index of the layer.
    */
    bool isConsensusFeatureVisible_(const Plot2DCanvas* canvas, const ConsensusFeature& cf, Size layer_index);

    const LayerDataConsensus* layer_; ///< the data to paint
  };

  /**
   @brief Painter2D for Identifications
  */
  class OPENMS_GUI_DLLAPI Painter2DIdent : public Painter2DBase
  {
  public:
    /// C'tor which remembers the layer to paint
    Painter2DIdent(const LayerDataIdent* parent);

    /// Implementation of base class
    void paint(QPainter*, Plot2DCanvas* canvas, int layer_index) override;

  protected:
    const LayerDataIdent* layer_; ///< the data to paint
  };
  
} // namespace OpenMS
