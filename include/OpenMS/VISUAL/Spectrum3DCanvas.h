// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public2
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
// $Maintainer: $
// $Authors: Cornelia Friedle $
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
    @brief Canvas for 3D-visualization of peak map data
		
		The Spectrum3DCanvas uses the helper class Spectrum3DOpenGLCanvas for the actual 3D rendering.
		Deriving Spectrum3DCanvas directly from QGLWidget is not possible due to the "Deadly Diamond" shape
		of inheritence.
		
		@image html Spectrum3DWidget.png
		
		@htmlinclude OpenMS_Spectrum3DCanvas.parameters
		
    @ingroup SpectrumWidgets
  */	
  class OPENMS_DLLAPI Spectrum3DCanvas
  	: public SpectrumCanvas
  {
    Q_OBJECT
    
    friend class Spectrum3DOpenGLCanvas;
    
	  public:
			
	    /// Constructor
	    Spectrum3DCanvas(const Param& preferences, QWidget* parent = 0);	
	    /// Destructor
	    virtual  ~Spectrum3DCanvas();
	    
	  	///Different shade modes
	    enum ShadeModes 
	    {
				SHADE_FLAT = 0,            
				SHADE_SMOOTH = 1         
	    };
	    
	    ///returns the Spectrum3DOpenGLcanvas     
	    Spectrum3DOpenGLCanvas* openglwidget();
	    
	    ///@name Remplemented Qt events
	    //@{
	    void resizeEvent(QResizeEvent* e);
			void contextMenuEvent(QContextMenuEvent* e);
	    //@}
	    /// Returns if the legend is shown
	    bool isLegendShown() const;
	    ///Shows/hides the legend
	    void showLegend(bool);
	    ///pointer to the SpectrumOpenGLCanvas implementation
	    Spectrum3DOpenGLCanvas* openglcanvas_;
			
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
  		
			/// Reacts on changed layer paramters
			void currentLayerParamtersChanged_();
  	
  	protected:
			
			// Docu in base class
	    bool finishAdding_();
			
  		// Reimplementation in order to update the OpenGL widget
      virtual void update_(const char* caller_name = 0);

			///whether the legend is shoen or not
			bool legend_shown_;

      //docu in base class
			virtual void translateLeft_();
			//docu in base class
			virtual void translateRight_();
			//docu in base class
			virtual void translateForward_();
			//docu in base class
			virtual void translateBackward_();
  };
  
} //namespace
#endif
