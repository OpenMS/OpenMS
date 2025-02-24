// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Cornelia Friedle $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QOpenGLWidget>
#include <QOpenGLFunctions_2_0>

// OpenMS
#include <OpenMS/DATASTRUCTURES/DRange.h>

namespace OpenMS
{
  class Plot3DCanvas;
  class LayerDataBase;

  /**
      @brief OpenGL Canvas for 3D-visualization of map data

      @note Do not use this class directly. Use Plot3DCanvas instead!

      @ingroup PlotWidgets
  */

  class OPENMS_GUI_DLLAPI Plot3DOpenGLCanvas :
    public QOpenGLWidget, 
    protected QOpenGLFunctions_2_0
  {
    Q_OBJECT

    friend class Plot3DCanvas;

public:

    /// Container for axis ticks
    typedef std::vector<std::vector<double> > AxisTickVector;

    /**
     @brief Constructor

     @param parent The parent widget
     @param canvas_3d The main 3d canvas
    */
    Plot3DOpenGLCanvas(QWidget * parent, Plot3DCanvas & canvas_3d);
    /**
        @brief Destructor

        Destroys the OpenGLWidget and all associated data.
    */
    ~Plot3DOpenGLCanvas() override;

    ///virtual function provided from QGLWidget
    void initializeGL() override;
    /// virtual function provided from QGLWidget
    void resizeGL(int w, int h) override;
    /// virtual function provided from QGLWidget
    void paintGL() override;

    /** @name Reimplemented QT events */
    //@{
    void mouseMoveEvent(QMouseEvent * e) override;
    void mouseReleaseEvent(QMouseEvent * e) override;
    void mousePressEvent(QMouseEvent * e) override;
    void focusOutEvent(QFocusEvent * e) override;
    //@}

    void setXLabel(const QString& l) { x_label_ = l; }
    void setYLabel(const QString& l) { y_label_ = l; }
    void setZLabel(const QString& l) { z_label_ = l; }
    
    /// updates the min and max values of the intensity
    void updateIntensityScale();
protected:
    /// helper function to project point to device space
    GLint project_(GLdouble objx, GLdouble objy, GLdouble objz, GLdouble * winx, GLdouble * winy); 
    /// helper function to transform point using matrix m (homogeneous coordinates)
    void transformPoint_(GLdouble out[4], const GLdouble m[16], const GLdouble in[4]);
    ///helper function to replicate old behaviour of QGLWidget
    void renderText_(double x, double y, double z, const QString & text);
    ///helper function to replicate old behaviour of QGLWidget
    void qglColor_(const QColor& color);
    ///helper function to replicate old behaviour of QGLWidget
    void qglClearColor_(const QColor& clearColor);
    /// Builds up a display list for the 3D view
    GLuint makeDataAsStick_();
    /// Builds up a display list for the axes
    GLuint makeAxes_();
    /// Builds up a display list for axis ticks
    GLuint makeAxesTicks_();
    /// Builds up a display list for the birds-eye view
    GLuint makeDataAsTopView_();
    /// Builds up a display list for the background
    GLuint makeGround_();
    /// Builds up a display list for grid lines
    GLuint makeGridLines_();
    /// Draws the axis texts
    void drawAxesLegend_();

    /// computes the dataset supposed to be drawn when a section has been selected in zoom mode
    void computeSelection_();

    /// calculates the zoom area , which is shown
    void dataToZoomArray_(double x_1, double y_1, double x_2, double y_2);

    /// returns the BB-rt-coordinate :  value --> BB-coordinates
    double scaledRT_(double rt);
    /// returns the rt-value : BB-coordinates  --> value
    double scaledInversRT_(double mz);
    /// returns the BB-mz-coordinate :  values --> BB-coordinates
    double scaledMZ_(double mz);
    ///  returns the mz-value : BB-coordinates  --> value
    double scaledInversMZ_(double mz);
    /// returns the BB-intensity -coordinate :  values --> BB-coordinates
    double scaledIntensity_(float intensity, Size layer_index);
    
    /// recalculates the dot gradient interpolation values.
    void recalculateDotGradient_(LayerDataBase& layer);
    ///calculate the ticks for the gridlines
    void calculateGridLines_();

    /// normalize the angle by "angle % 360*16"
    void normalizeAngle(int* angle);
    // set translation vector to 0
    void resetTranslation();

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

    /// reference to Plot3DCanvas
    Plot3DCanvas & canvas_3d_;

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

    QPainter* painter_ = nullptr;

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

    QString x_label_;
    QString y_label_;
    QString z_label_;

protected slots:
    /// Slot that reacts on action mode changes
    void actionModeChange();
  };
}
