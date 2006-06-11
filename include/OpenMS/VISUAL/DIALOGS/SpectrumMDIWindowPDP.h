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
// $Id: SpectrumMDIWindowPDP.h,v 1.11 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUMMDIWINDOWPDP_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUMMDIWINDOWPDP_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialogPage.h>
#include <OpenMS/VISUAL/SpectrumMDIWindow.h>

class QComboBox;
class QLineEdit;
class QRadioButton;
class QCheckBox;
class QSpinBox;

namespace OpenMS
{
	class MultiGradientSelector;
	class ColorSelector;
	
	namespace Internal
	{
		///Preferences dialog page of a SpectrumMDIWindow (internal use only)	
		class SpectrumMDIWindowPDP: public PreferencesDialogPage
		{
			Q_OBJECT
			
		public:
			SpectrumMDIWindowPDP( SpectrumMDIWindow* manager, QWidget* parent = 0, const char* name = "SpectrumMDIWindowPDP", WFlags f = 0);
			virtual ~SpectrumMDIWindowPDP();
			virtual void load();
			virtual void save();
		
		public slots:
			void browseDefaultPath_();
		
		protected:
			QLineEdit* main_default_path_;
			QSpinBox* recent_files_;
			QRadioButton* map_view_2d_;
			QRadioButton* map_view_3d_;
			
			//db
			QLineEdit* db_host_;
			QLineEdit* db_port_;
			QLineEdit* db_name_;
			QLineEdit* db_login_;
			
			//1d
			ColorSelector* peak_color_;
			ColorSelector* high_color_;
			ColorSelector* icon_color_;
			ColorSelector* back_color_1D_;
			ColorSelector* back_color_2D_;
			QComboBox* axis_mapping_2d_;
			QComboBox* axis_mapping_;
			QComboBox* x_axis_orientation_;
			QComboBox* y_axis_orientation_;
			QCheckBox* log_check_box_;
			QCheckBox* rel_check_box_;
			QCheckBox* x_show_legend_;
			QCheckBox* y_show_legend_;
			
			//2d
			QSpinBox* marching_squares_steps_;
			QSpinBox* dot_interpolation_steps_;
			QSpinBox* surface_interpolation_steps_;
			QRadioButton* dot_mode_black_;
			QRadioButton* dot_mode_gradient_;
			MultiGradientSelector* dot_gradient_;
			MultiGradientSelector* surface_gradient_;
			QComboBox* x_axis_orientation_2d_;
			QComboBox* y_axis_orientation_2d_;
			
			//3d
		
			QRadioButton* dot_mode_black_3d_;
			QRadioButton* dot_mode_gradient_3d_;
			ColorSelector* back_color_3d_; 
			MultiGradientSelector* dot_gradient_3d_;
			QSpinBox* dot_interpolation_steps_3d_;
			QRadioButton* shade_mode_flat_3d_;
			QRadioButton* shade_mode_smooth_3d_;
			QRadioButton* intensity_mode_lin_3d_;
			QRadioButton* intensity_mode_log_3d_;
			ColorSelector* axes_color_3d_;
		};
	
	} //namespace Internal

	
} //namespace OpenMS

#endif //OPENMS_VISUAL_DIALOGS_SPECTRUMMDIWINDOWPDP_H

