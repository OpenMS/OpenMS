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
   
     @todo Fix taking of images (Cornelia)
     
     @todo Axis labels (Cornelia)
     
     @todo Add translation

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
    
    // Docu in base class
    virtual PreferencesDialogPage* createPreferences(QWidget* parent);
    
    ///returns the Spectrum3DOpenGLcanvas     
    Spectrum3DOpenGLCanvas* openglwidget();
    
    // Docu in base class
    virtual void invalidate_();
		// Docu in base class
    SignedInt finishAdding(float low_intensity_cutoff = 0);
	  // Docu in base class
    virtual void actionModeChange_();
	  // Docu in base class
    virtual void setMainPreferences(const Param& prefs);
		// Docu in base class
		virtual void intensityModeChange_();
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
