// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM3DCANVAS_H

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>

class QPainter;
class QGLWidget;
class QResizeEvent;

namespace OpenMS
{
  class Spectrum3DOpenGLCanvas;
  
  /**
    @brief Canvas for 3D-visualization of map data
		
		@ref Spectrum3DCanvas_Parameters are explained on a separate page.
		
    @ingroup SpectrumWidgets
  */	
  class Spectrum3DCanvas
  	: public SpectrumCanvas
  {
    Q_OBJECT
    
    friend class Spectrum3DOpenGLCanvas;
    
	  public:
	    /// Constructor
	    Spectrum3DCanvas(const Param& preferences, QWidget* parent = 0);	
	    /// Destructor
	    virtual  ~Spectrum3DCanvas();
	    
	    /**	@name Type definitions */
	    //@{
	     ///different reduktion modes: no reduction, sumreduktion or maxreduktion
 			enum DataModes 
 	    {
 				REDUCTION_OFF = 0,          
 				REDUCTION_MAX = 1,
 				REDUCTION_SUM = 2
 	    };
	  	///Different shade modes
	    enum ShadeModes 
	    {
				SHADE_FLAT = 0,            
				SHADE_SMOOTH = 1         
	    };
	    //@}
	    
	    ///returns the Spectrum3DOpenGLcanvas     
	    Spectrum3DOpenGLCanvas* openglwidget();

			// Docu in base class
	    Int finishAdding();
	    
	    ///QT resizeEvent
	    void resizeEvent(QResizeEvent * e);
	    /// Returns if the legend is shown
	    bool isLegendShown() const;
	    ///Shows/hides the legend
	    void showLegend(bool);
	    ///pointer to the SpectrumOpenGLCanvas implementation
	    Spectrum3DOpenGLCanvas* openglcanvas_;
			///makes the reduced datasets using sumreducer or maxreducer
	    void makeReducedDataSet();
			
			// docu in base class
			virtual void showCurrentLayerPreferences();
			
		public slots:
	    // Docu in base class
	    void activateLayer(int layer_index);
	    // Docu in base class
	    void removeLayer(int layer_index);
  	
  	protected:
  		// Docu in base class
			virtual void intensityModeChange_();
			/**
				@brief Sets the visible area
				
				Changes the visible area, adjustes the zoom stack and notifies interested clients about the change. 
				If parts of the area are outside of the data area, the new area will be adjusted.
				
				@param new_area The new visible area.
				@param add_to_stack If the new area is to add to the zoom_stack_
			*/
			virtual void changeVisibleArea_(const AreaType& new_area, bool add_to_stack = false);
			
			// Docu in base class
			virtual void currentLayerParamtersChanged_();
			
  		// Reimplementation in order to update the OpenGL widget
  		virtual void update_(const char* caller_name);

			/// area of the layer
			double area_;
			///whether the legend is shoen or not
			bool legend_shown_;
  };
  
} //namespace
#endif
