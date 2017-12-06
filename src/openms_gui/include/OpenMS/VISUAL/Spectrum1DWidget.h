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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM1DWIDGET_H
#define OPENMS_VISUAL_SPECTRUM1DWIDGET_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

// STL
#include <vector>

// OpenMS
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

class QAction;
class QSpacerItem;

namespace OpenMS
{
  class Spectrum1DCanvas;

  /**
      @brief Widget for visualization of several spectra

      The widget is composed of a scroll bar, an AxisWidget and a Spectrum1DCanvas as central widget.

      @image html Spectrum1DWidget.png

      The example image shows %Spectrum1DWidget displaying a raw data layer and a peak data layer.

      @ingroup SpectrumWidgets
  */
  class OPENMS_GUI_DLLAPI Spectrum1DWidget :
    public SpectrumWidget
  {
    Q_OBJECT

public:
    /// Default constructor
    Spectrum1DWidget(const Param & preferences, QWidget * parent = nullptr);
    ///Destructor
    ~Spectrum1DWidget() override;

    /// This method is overwritten to make the class specific members accessible
    inline Spectrum1DCanvas * canvas()
    {
      return static_cast<Spectrum1DCanvas *>(canvas_);
    }

    // Docu in base class
    void hideAxes() override;

    // Docu in base class
    void showLegend(bool show) override;

    /// Switches to mirror view, displays another y-axis for the second spectrum
    void toggleMirrorView(bool mirror);

    /// Performs an alignment of the layers with @p layer_index_1 and @p layer_index_2
    void performAlignment(Size layer_index_1, Size layer_index_2, const Param & param);

    /// Resets the alignment
    void resetAlignment();

    // Docu in base class
    void saveAsImage() override;

    // Docu in base class
    virtual void renderForImage(QPainter& painter);

signals:
    /// Is emitted whenever the visible area changes.
    void visibleAreaChanged(double, double);

    /// Requests to display the whole spectrum in 2D view
    void showCurrentPeaksAs2D();

    /// Requests to display the whole spectrum in 3D view
    void showCurrentPeaksAs3D();

public slots:
    // Docu in base class
    void showGoToDialog() override;

protected:
    // Docu in base class
    Math::Histogram<> createIntensityDistribution_() const override;
    // Docu in base class
    Math::Histogram<> createMetaDistribution_(const String & name) const override;
    // Docu in base class
    void recalculateAxes_() override;

    /// The second y-axis for the mirror view
    AxisWidget * flipped_y_axis_;

    /// Spacer between the two y-axes in mirror mode (needed when visualizing an alignment)
    QSpacerItem * spacer_;

  };
} // namespace OpenMS

#endif
