
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
// $Maintainer:Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM3DCANVASPDP_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM3DCANVASPDP_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialogPage.h>
#include <qradiobutton.h>
class QSpinBox;

namespace OpenMS
{
  class ColorSelector;
  class Spectrum3DCanvas;
  class MultiGradientSelector;
  
  namespace Internal
  {
    
    ///Preferences dialog page of a Spectrum3DCanvas (internal use only)	
    class Spectrum3DCanvasPDP: public PreferencesDialogPage
    {
      Q_OBJECT
      
    public:
      Spectrum3DCanvasPDP( Spectrum3DCanvas* manager, QWidget* parent = 0, const char* name = "Spectrum3DWidgetPDP", WFlags f = 0);
      virtual ~Spectrum3DCanvasPDP();
      virtual void load();
      virtual void save();
      
    protected:
      
      
      QRadioButton* dot_mode_black_;
      QRadioButton* dot_mode_gradient_;
      ColorSelector* background_color_;
      MultiGradientSelector* dot_gradient_;
      QSpinBox* dot_interpolation_steps_;
      QRadioButton* shade_mode_flat_;
      QRadioButton* shade_mode_smooth_;
      
      QRadioButton* reduction_on_max_;
      QRadioButton* reduction_on_sum_;
      QRadioButton* reduction_off_;
      QSpinBox* reduction_ratio_max_;
      QSpinBox* reduction_ratio_sum_;
    
      ColorSelector* axes_color_;
      QSpinBox* dot_line_width_;
    };
    
  } //namespace Internal
  
} //namespace OpenMS

#endif 
