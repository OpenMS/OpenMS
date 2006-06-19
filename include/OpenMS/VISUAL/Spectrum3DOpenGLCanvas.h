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
// $Id: Spectrum3DOpenGLCanvas.h,v 1.27 2006/06/09 11:47:50 cfriedle Exp $
// $Author: cfriedle $
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DOPENGLCANVAS_H
#define OPENMS_VISUAL_SPECTRUM3DOPENGLCANVAS_H

// QT
#include <qpoint.h>
#include <qgl.h>

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
		typedef DSpectrum< 1, OpenMS::DPeakArrayNonPolymorphic<1, PeakT > > 	BaseSpectrum;
		//@}



 public:
    /**
     @brief Constructor
     @param parent The parent widget
     @param name The widget's name
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
		
		/// method to make the font
		void	makeFont();
		/// method to make the font		
		GLuint fontOffset;
		/// method to make the font
		void paintAxesScale();

		///Mouse-and Wheel- Events
    void mouseMoveEvent(QMouseEvent *e);
		void mouseReleaseEvent(QMouseEvent *e);
    void mousePressEvent(QMouseEvent *e);
		void wheelEvent(QWheelEvent * e);
		/**
			 @brief Sets the 3D stick gradient.
			 the dot mode is set to	DOT_GRADIENT.
			 @param gradient A string containing the gradient description.
		*/
		void setDotGradient(const std::string& gradient);
		///



		void updateMinMaxValues();
		void updateIntensityScale();
		void dataToZoomArray(double x_1, double y_1, double x_2, double y_2);
		///

		double scaledRT(double rt);
		double scaledZoomRT(double rt);
		double scaledInversRT(double mz);
	  double scaledMZ(double mz);
		double scaledZoomMZ(double mz);
		double scaledInversMZ(double mz);
    double scaledIntensity(double intensity);
		double scaledZoomIntensity(double intensity);
		double scaledLogIntensity(double intensity);
		/// recalculates the dot gradient inerpolation values.
		void recalculateDotGradient_(UnsignedInt);
		/// recalculates the logarthmic dot gradient inerpolation values.
		void recalculateDotGradientLog_(UnsignedInt i);
		
		///calculate the ticks for the gridlines
		void calculateGridLines_();
		
		
		//views from toolbar
		void setBackView();
		void setTopView();
		void setResetZoomView();
		void setZoomView();
		void setSelectView();
		bool getShowSelect();		
		bool getShowZoom();
		void setIntensityScale(bool);

    ///
    int xRotation() const { return xRot_; }
    int yRotation() const { return yRot_; }
    int zRotation() const { return zRot_; }
    void normalizeAngle(int *angle);
		
		ViewMode view_mode_;
		//displaylisten
    GLuint stickdata_;
		GLuint coord_;
		GLuint axeslabel_;
		GLuint zoomselection_;
		GLuint zoomdata_;
		GLuint ground_;

		//preferences
		MultiGradient gradient_;

		// reference to Spectrum3DCanvas
		Spectrum3DCanvas& canvas_3d_;
  
		/// member variables for the rotation
    int xRot_,xRot_old_;
    int yRot_,yRot_old_;
    int zRot_,zRot_old_;

		/// member variables fot the zoom-modus
    QPoint lastMousePos_,firstMousePos_;    

		///member vairables for the BB ans the resize event
    double corner_, zoom_ ,near_, far_;
		float width_;
    float heigth_;

		bool second_paint_;
		bool topview_,showbackview_;
		bool zoom_mode_;
		bool show_zoom_selection_;		
		bool grid_exists_;
		bool intensity_scale_;
		
		DRange<3> overall_values_;
		DRange<1> int_scale_;
		///member gridvectors which contains the data for the ticks
		GridVector grid_mz_,grid_rt_, grid_intensity_,grid_intensity_log_;

		double x_1_,x_2_,y_1_,y_2_;

public slots:
    void setRotationX(int angle);
    void setRotationY(int angle);
    void setRotationZ(int angle);
		void setZoomFactor(double zoom);
signals:
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void zRotationChanged(int angle);
		void zoomFactorChanged(double zoom);

		void rightButton(QPoint pos);
	};
  
  
}
#endif
