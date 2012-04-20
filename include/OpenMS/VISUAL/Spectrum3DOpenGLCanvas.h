// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DOPENGLCANVAS_H
#define OPENMS_VISUAL_SPECTRUM3DOPENGLCANVAS_H

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
    virtual ~Spectrum3DOpenGLCanvas();
    ///virtual function provided from QGLWidget
    void initializeGL();
    /// virtual function provided from QGLWidget
    void resizeGL(int w, int h);
    /// virtual function provided from QGLWidget
    void paintGL();
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
    void mouseMoveEvent(QMouseEvent * e);
    void mouseReleaseEvent(QMouseEvent * e);
    void mousePressEvent(QMouseEvent * e);
    void focusOutEvent(QFocusEvent * e);
    //@}

    /// computes the dataset supposed to be drawn when a section has been selected in zoom mode
    void computeSelection();

    /// updates the min and max values of the intensity
    void updateIntensityScale();

    /// calcualtes the zoom area , which is shown
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
    double scaledIntensity(Real intensity, Size layer_index);

    /// recalculates the dot gradient inerpolation values.
    void recalculateDotGradient_(Size layer);
    ///calculate the ticks for the gridlines
    void calculateGridLines_();

    /// return width
    float width() const { return width_; }
    float height() const { return heigth_; }

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



    /// member variables fot the zoom-modus
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
    float heigth_;
    /// object which contains the min and max values of mz, rt and intensity
    DRange<3> overall_values_;
    ///object wich contains the values of the current min and max intensity
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
