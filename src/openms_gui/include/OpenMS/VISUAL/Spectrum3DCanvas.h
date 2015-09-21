// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM3DCANVAS_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>

class QPainter;
class QGLWidget;
class QResizeEvent;

namespace OpenMS
{
  class Spectrum3DOpenGLCanvas;

  /**
    @brief Canvas for 3D-visualization of peak map data

        The Spectrum3DCanvas uses the helper class Spectrum3DOpenGLCanvas for the actual 3D rendering.
        Deriving Spectrum3DCanvas directly from QGLWidget is not possible due to the "Deadly Diamond" shape
        of inheritance.

        @image html Spectrum3DWidget.png

        @htmlinclude OpenMS_Spectrum3DCanvas.parameters

    @ingroup SpectrumWidgets
  */
  class OPENMS_GUI_DLLAPI Spectrum3DCanvas :
    public SpectrumCanvas
  {
    Q_OBJECT

    friend class Spectrum3DOpenGLCanvas;

public:

    /// Constructor
    Spectrum3DCanvas(const Param & preferences, QWidget * parent = 0);
    /// Destructor
    virtual  ~Spectrum3DCanvas();

    ///Different shade modes
    enum ShadeModes
    {
      SHADE_FLAT = 0,
      SHADE_SMOOTH = 1
    };

    ///returns the Spectrum3DOpenGLcanvas
    Spectrum3DOpenGLCanvas * openglwidget();

    ///@name Reimplemented Qt events
    //@{
    void resizeEvent(QResizeEvent * e);
    void contextMenuEvent(QContextMenuEvent * e);
    //@}
    /// Returns if the legend is shown
    bool isLegendShown() const;
    ///Shows/hides the legend
    void showLegend(bool);
    ///pointer to the SpectrumOpenGLCanvas implementation
    Spectrum3DOpenGLCanvas * openglcanvas_;

    // docu in base class
    virtual void showCurrentLayerPreferences();

    // Docu in base class
    virtual void saveCurrentLayer(bool visible);

signals:
    /// Requests to display all spectra in 2D plot
    void showCurrentPeaksAs2D();

public slots:

    // Docu in base class
    void activateLayer(Size layer_index);
    // Docu in base class
    void removeLayer(Size layer_index);
    //docu in base class
    virtual void updateLayer(Size i);
protected slots:

    /// Reacts on changed layer parameters
    void currentLayerParamtersChanged_();

protected:

    // Docu in base class
    bool finishAdding_();

    // Reimplementation in order to update the OpenGL widget
    virtual void update_(const char * caller_name = 0);

    ///whether the legend is shown or not
    bool legend_shown_;

    //docu in base class
    virtual void translateLeft_(Qt::KeyboardModifiers m);
    //docu in base class
    virtual void translateRight_(Qt::KeyboardModifiers m);
    //docu in base class
    virtual void translateForward_();
    //docu in base class
    virtual void translateBackward_();
  };

} //namespace
#endif
