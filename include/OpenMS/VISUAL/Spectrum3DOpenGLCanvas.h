// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DOPENGLCANVAS_H
#define OPENMS_VISUAL_SPECTRUM3DOPENGLCANVAS_H

// QT
#include <qpoint.h>
#include <qpixmap.h>
#include <qgl.h>

#include<qpainter.h>
// STL
#include <vector>

// OpenMS
#include <OpenMS/config.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include<OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
namespace OpenMS
{

	class Spectrum3DCanvas;
	class Param;
	
	/**
		@brief OpenGL Canvas for 3D-visualization of map data
		
		@ingroup spectrum_widgets
	*/
	
  class Spectrum3DOpenGLCanvas:public QGLWidget
  {
    Q_OBJECT
		
		friend class Spectrum3DCanvas;
		
		/**	@name Type definitions */
		//@{
		enum ViewMode
		{
			VIEW_SELECT,
			VIEW_TOP,
			VIEW_ZOOM
		};
		///
		typedef std::vector<std::vector<double> > GridVector;
		///
		typedef DimensionDescription < DimensionDescriptionTagLCMS > DimDesc;
		///
		enum DimensionId { MZ = DimDesc::MZ, RT = DimDesc::RT };
		
		typedef DPeak<1> PeakT;
		typedef DSpectrum< 1, DPeakArrayNonPolymorphic<1, PeakT > > 	BaseSpectrum;
		//@}

 public:
    /**
     @brief Constructor
     
     @param parent The parent widget
     @param name The widget's name
     @param canvas_3d The main 3d canvas
    */
		Spectrum3DOpenGLCanvas(QWidget *parent, const char* name, Spectrum3DCanvas & canvas_3d);
		/**
			@brief Destructor
		
			Destroys the OpenGLWidget and all associated data.
		*/		
    virtual ~Spectrum3DOpenGLCanvas();
		/**
			 @brief virtual function provided from QGLWidget	
		*/
    void initializeGL();
		/**
			 @brief virtual function provided from QGLWidget	
		*/ 
		void resizeGL(int w,int h);
    /**
			 @brief virtual function provided from QGLWidget	
		*/
		void paintGL();
		/**
			 @brief Display list for the sticks which display the 3-dimensional data
		*/
    virtual GLuint makeDataAsStick();
	  /**
		  @brief Display list for the coordinates
	  */ 
		virtual GLuint makeCoordinates();
		/**
			 @brief Display list for the different axes label
		*/
		virtual GLuint makeAxesLabel();
	  /**
		 	 @brief Display list for the sticks 
		*/
		virtual GLuint makeDataAsStickLog();
		/**
			 @brief Display list for the top view
		*/
		virtual GLuint makeDataAsTopView();
		/**
			 @brief Display list 
		*/
		virtual GLuint makeZoomSelection();
		/**
		 * @brief displaylist
		 */	
		virtual GLuint makeGround();
		/**
		 * @brief displaylist
		 */	
		virtual GLuint makeGridLines();
		/// method to make the font
		virtual GLuint makeLegend();

		///reimplementation if the mouseMoveEvent
    void mouseMoveEvent(QMouseEvent *e);
		///reimplementation if the mousePressEvent
		void mouseReleaseEvent(QMouseEvent *e);
		///reimplementation if the mouseRealeaseEvent
    void mousePressEvent(QMouseEvent *e);
		///reimplementation if the wheelEvent
		void wheelEvent(QWheelEvent * e);
		void keyPressEvent(QKeyEvent * e);
		/**
			 @brief Sets the 3D stick gradient.
			 the dot mode is set to	DOT_GRADIENT.
			 @param gradient A string containing the gradient description.
		*/
		void setDotGradient(const std::string& gradient);

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
		double scaledIntensity(double intensity,int dataset);

		/// recalculates the dot gradient inerpolation values.
		void recalculateDotGradient_();
		///calculate the ticks for the gridlines
		void calculateGridLines_();
	
    /// return xRot_
    int xRotation() const { return xrot_; }
    /// return yRot_
		int yRotation() const { return yrot_; }
    /// return zRot_
		int zRotation() const { return zrot_; }
    /// normalize the angel
		void normalizeAngle(int *angle);
		///
		void resetAngels();
		///
		void resetTranslation();
		///
		ViewMode view_mode_;
		/// displaylist
    GLuint stickdata_;
		/// displaylist
		GLuint coord_;
		/// displaylist
		GLuint axeslabel_;
		///displaylist
		GLuint gridlines_;
		/// displaylist
		GLuint zoomselection_;
		/// displaylist
		GLuint zoomdata_;
		/// displaylist
		GLuint ground_;
		/// displaylist
		GLuint axeslegend_;
		//preferences
		MultiGradient gradient_;

		// reference to Spectrum3DCanvas
		Spectrum3DCanvas& canvas_3d_;
  
		/// member x-variables for the rotation
    int xrot_;
		/// member y-variables for the rotation
		int yrot_;
		/// member z-variables for the rotation
    int zrot_;

		/// member variables fot the zoom-modus
    QPoint lastMousePos_,firstMousePos_;    

		///member variable for the x and y axis of the BB 
    double corner_;
		/// member variable for the zoom- Modus
		double zoom_ ;
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
		GridVector grid_mz_;
		///member gridvectors which contains the data for the rt-axis-ticks
		GridVector grid_rt_;
		///member gridvectors which contains the data for the intensity-axis-ticks
		GridVector grid_intensity_;
		///member gridvectors which contains the data for the log-intensity-axis-ticks
		GridVector grid_intensity_log_;
		/// if it is true the zoom_selection is shown
		bool show_zoom_selection_;	
		///if the ControlKey is pressed ->true
		bool translation_on_;	
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
	// 	QPixmap *pix_;
		QPainter paint_;
		bool zoom_paint_;
public slots:
    /// first normalize the angel and then set xRot_ 
    void setRotationX(int angle);
		/// first normalize the angel and then set yRot_ 
    void setRotationY(int angle);
		/// first normalize the angel and then set zRot_ 
    void setRotationZ(int angle);
		/// set the member variable zoom_ and calls initializeGL and updateGL
		void setZoomFactor(double zoom);
signals:
		/// signal which is emmited in the mouseReleaseEvent 
		void rightButton(QPoint pos);
	};
  
  
}
#endif
