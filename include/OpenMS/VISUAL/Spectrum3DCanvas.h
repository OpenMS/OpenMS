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
// $Id: Spectrum3DCanvas.h,v 1.26 2006/06/09 11:53:40 cfriedle Exp $
// $Author: cfriedle $
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM3DCANVAS_H

// QT
#include <qpoint.h>

// STL

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/config.h>
#include <OpenMS/KERNEL/DimensionDescription.h>

class QPainter;
class QGLWidget;

namespace OpenMS
{
  class Spectrum3DOpenGLCanvas;
  
  /**
     @brief Canvas for 3D-visualization of map data
   
     @todo Fix intensity distribution (Cornelia)
     
     @todo Fix taking of images (Cornelia)
     
     @todo Replace IntScale by SpectrumCanvas::IntensityModes (Cornelia)
     
     @todo Axis labels (Cornelia)
     
     @todo Add translation
     
     @todo /share/data/MS/Examples/mzData/MS_small.mzData crashes when opening
     
     @todo Make axes visible in zoom/top view
     
     @todo mini.dta2d
     
     @ingroup spectrum_widgets
  */	
  class Spectrum3DCanvas : public SpectrumCanvas
  {
    Q_OBJECT
    
    friend class Spectrum3DOpenGLCanvas;
    
  public:
    /// Constructor
    Spectrum3DCanvas(QWidget* parent = 0, const char* name = "Spectrum3DCanvas", WFlags f = 0);	
    /// Destructor
    virtual  ~Spectrum3DCanvas();
    
    /**	@name Type definitions */
    //@{
	
    ///
    typedef DimensionDescription < DimensionDescriptionTagLCMS > DimDesc;
    ///
    enum DimensionId { MZ = DimDesc::MZ, RT = DimDesc::RT };	
    ///
    enum DotModes 
      {
			DOT_BLACK = 0,            ///< use black only
			DOT_GRADIENT = 1          ///< use gradient
      };
    ///
    enum ShadeModes 
      {
			SHADE_FLAT = 0,            
			SHADE_SMOOTH = 1         
      };
    ///
    enum IntScale 
      {
			INT_LINEAR = 0,            
			INT_LOG = 1         
      };
    //@}
    
    /**
       @brief Creates a preferences dialog page
	       
       
       @param parent the parent widget for the dialog page
    */
    virtual PreferencesDialogPage* createPreferences(QWidget* parent);
    
    ///returns the Spectrum3DOpenGLcanvas     
    Spectrum3DOpenGLCanvas* openglwidget();
    
    // Docu in base class
    virtual void invalidate_();
		// Docu in base class
    SignedInt finishAdding();
	  // Docu in base class
    virtual void actionModeChange_();
	  // Docu in base class
    virtual void setMainPreferences(const Param& prefs);
        
    ///preferences
    SignedInt getDotMode();
    void setDotGradient(const std::string& gradient);
    SignedInt getShadeMode();
    UnsignedInt getDotInterpolationSteps();
    
    //resizeEvent
    void resizeEvent(QResizeEvent * e);
    void connectMouseEvents();
    
    Spectrum3DOpenGLCanvas* openglcanvas_;
    int current_zoom_;
    
public slots:
    ///shows the contextmenu at position p
    void showContextMenu(QPoint p);
    // Docu in base class
    void activateDataSet(int data_set);
    // Docu in base class
    void removeDataSet(int data_set);
  };
  
} //namespace
#endif
