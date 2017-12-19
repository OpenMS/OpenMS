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
// $Maintainer: Timo Sachsenberg $
// $Authors: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DOPENGLCANVAS_H
#define OPENMS_VISUAL_SPECTRUM3DOPENGLCANVAS_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtOpenGL/QGLWidget>

// OpenMS
#include <OpenMS/DATASTRUCTURES/DRange.h>

namespace OpenMS
{
  class Spectrum3DCanvas;

  /**
      @brief OpenGL Canvas for 3D-visualization of map data

      @note Do not use this class directly. Use Spectrum3DCanvas instead!

      @ingroup SpectrumWidgets
  */

  class OPENMS_GUI_DLLAPI Spectrum3DOpenGLCanvas :
    public QGLWidget
  {
    Q_OBJECT

    friend class Spectrum3DCanvas;

public:

    /// Container for axis ticks
    typedef std::vector<std::vector<double> > AxisTickVector;

    /**
     @brief Constructor

     @param parent The parent widget
     @param canvas_3d The main 3d canvas
    */
    Spectrum3DOpenGLCanvas(QWidget * parent, Spectrum3DCanvas & canvas_3d);
    /**
        @brief Destructor

        Destroys the OpenGLWidget and all associated data.
    */
    ~Spectrum3DOpenGLCanvas() override;
    ///virtual function provided from QGLWidget
    void initializeGL() override;
    /// virtual function provided from QGLWidget
    void resizeGL(int w, int h) override;
    /// virtual function provided from QGLWidget
    void paintGL() override;
    /// Builds up a display list for the 3D view
    GLuint makeDataAsStick();
    /// Builds up a display list for the axes
    GLuint makeAxes();
    /// Builds up a display list for axis ticks
    GLuint makeAxesTicks();
    /// Builds up a display list for the birds-eye view
    GLuint makeDataAsTopView();
    /// Builds up a display list for the background
    GLuint makeGround();
    /// Builds up a display list for grid lines
    GLuint makeGridLines();
    /// Draws the axis texts (since Qt 4.3 these cannot be put into display lists anymore...)
    void drawAxesLegend();

    /** @name Reimplemented QT events */
    //@{
    void mouseMoveEvent(QMouseEvent * e) override;
    void mouseReleaseEvent(QMouseEvent * e) override;
    void mousePressEvent(QMouseEvent * e) override;
    void focusOutEvent(QFocusEvent * e) override;
    //@}

    /// computes the dataset supposed to be drawn when a section has been selected in zoom mode
    void computeSelection();

    /// updates the min and max values of the intensity
    void updateIntensityScale();

    /// calculates the zoom area , which is shown
    void dataToZoomArray(double x_1, double y_1, double x_2, double y_2);

    /// returns the BB-rt-coordinate :  value --> BB-coordinates
    double scaledRT(double rt);
    /// returns the rt-value : BB-coordinates  --> value
    double scaledInversRT(double mz);
    /// returns the BB-mz-coordinate :  values --> BB-coordinates
    double scaledMZ(double mz);
    ///  returns the mz-value : BB-coordinates  --> value
    double scaledInversMZ(double mz);
    /// returns the BB-intensity -coordinate :  values --> BB-coordinates
    double scaledIntensity(float intensity, Size layer_index);

    /// recalculates the dot gradient interpolation values.
    void recalculateDotGradient_(Size layer);
    ///calculate the ticks for the gridlines
    void calculateGridLines_();

    /// return width
    float width() const { return width_; }
    float height() const { return height_; }

    /// return xRot_
    int xRotation() const { return xrot_; }
    /// return yRot_
    int yRotation() const { return yrot_; }
    /// return zRot_
    int zRotation() const { return zrot_; }
    /// normalize the angel
    void normalizeAngle(int * angle);
    //document me
    void setAngels(int xrot, int yrot, int zrot);
    //document me
    void resetTranslation();
    //document me
    void timeMessure();

    /// stores the original rotation and zoom factor (e.g. before changing into zoom mode)
    void storeRotationAndZoom();
    /// restores the original rotation and zoom factor (e.g. before changing into zoom mode)
    void restoreRotationAndZoom();

    /** @name Different OpenGL display lists */
    //@{
    GLuint stickdata_;
    GLuint axes_;
    GLuint axes_ticks_;
    GLuint gridlines_;
    GLuint ground_;
    //@}

    /// reference to Spectrum3DCanvas
    Spectrum3DCanvas & canvas_3d_;

    /// member x-variables for the rotation
    int xrot_;
    /// member y-variables for the rotation
    int yrot_;
    /// member z-variables for the rotation
    int zrot_;

    /// member x-variable that stores the original angle during zoom mode
    int xrot_tmp_;
    /// member y-variable that stores the original angle during zoom mode
    int yrot_tmp_;
    /// member z-variable that stores the original angle during zoom mode
    int zrot_tmp_;



    /// member variables for the zoom-mode
    QPoint mouse_move_end_, mouse_move_begin_;

    ///member variable for the x and y axis of the BB
    double corner_;
    /// member variable for the zoom mode
    double zoom_;
    /// member variable that stores original zoom factor during zoom mode
    double zoom_tmp_;

    /// member variable for the z- axis of the BB
    double near_;
    /// member variable for the z- axis of the BB
    double far_;
    /// the width of the viewport
    float width_;
    /// the height of the viewport
    float height_;
    /// object which contains the min and max values of mz, rt and intensity
    DRange<3> overall_values_;
    ///object which contains the values of the current min and max intensity
    DRange<1> int_scale_;
    ///member gridvectors which contains the data for the mz-axis-ticks
    AxisTickVector grid_mz_;
    ///member gridvectors which contains the data for the rt-axis-ticks
    AxisTickVector grid_rt_;
    ///member gridvectors which contains the data for the intensity-axis-ticks
    AxisTickVector grid_intensity_;
    /// x1 coordinate of the zoomselection
    double x_1_;
    /// x2 coordinate of the zoomselection
    double x_2_;
    /// y1 coordinate of the zoomselection
    double y_1_;
    /// y2 coordinate of the zoomselection
    double y_2_;
    /// x- translation
    double trans_x_;
    /// y_translation
    double trans_y_;

protected slots:
    /// Slot that reacts on action mode changes
    void actionModeChange();
  };
}
#endif
