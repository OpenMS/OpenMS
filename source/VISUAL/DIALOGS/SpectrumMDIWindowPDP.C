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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/VISUAL/DIALOGS/SpectrumMDIWindowPDP.h>
#include <OpenMS/VISUAL/SpectrumMDIWindow.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/PreferencesManager.h>
#include <OpenMS/VISUAL/ColorSelector.h>
// Qt
#include <qlayout.h>
#include <qlineedit.h>
#include <qradiobutton.h>
#include <qtabwidget.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qcheckbox.h>
#include <qvbuttongroup.h>
#include <qhbuttongroup.h>
#include <qfiledialog.h>
#include <qcombobox.h>
#include <qpushbutton.h>
#include <qspinbox.h>

using namespace std;

namespace OpenMS
{

	namespace Internal
	{


		SpectrumMDIWindowPDP::SpectrumMDIWindowPDP( SpectrumMDIWindow* manager, QWidget* parent, const char* name, WFlags f)
			:PreferencesDialogPage(manager,parent,name,f)
		{
			help_ = "This is the preferences dialog of the main window!"
							"<br>";
		
			QGridLayout* grid;
			QLabel* label;
			QWidget* background;
		
			//tab widget
			grid = new QGridLayout(this);
			QTabWidget* tab = new QTabWidget(this);
			grid->addWidget(tab,0,0);
			
			//-----------General Tab-----------
			background = new QWidget(this);
			grid = new QGridLayout(background,2,2);
			grid->setMargin(6);
			grid->setSpacing(4);
			
			//default path
			label = new QLabel("Default path:",background);
			grid->addWidget(label,0,0,AlignLeft);
			main_default_path_ = new QLineEdit(background);
			grid->addWidget(main_default_path_,0,1,AlignLeft);
			main_default_path_->setMinimumWidth(fontMetrics().width('W') * 25);
			QPushButton* tmp = new QPushButton("Browse",background);
			grid->addWidget(tmp,0,2,AlignLeft);
			connect(tmp,SIGNAL(clicked()),this,SLOT(browseDefaultPath_()));
			//recent files
			label = new QLabel("Number of recent files:",background);
			grid->addWidget(label,1,0,AlignLeft);
			recent_files_ = new QSpinBox(1,99,1,background);
			grid->addWidget(recent_files_,1,1,AlignLeft);
			//default map view
			label = new QLabel("Default map visualization:",background);
			grid->addWidget(label,2,0,AlignTop);
			default_map_view_ = new QComboBox(false, background, "read-only combobox");
			default_map_view_->insertItem("2D");
			default_map_view_->insertItem("3D");
			grid->addWidget(default_map_view_,2,1,AlignLeft);
			//legend
			label = new QLabel("Axis legend: ",background);
			grid->addWidget(label,3,0,AlignTop);
			show_legend_ = new QComboBox(false, background, "read-only combobox");
			show_legend_->insertItem("Show");
			show_legend_->insertItem("Hide");
			grid->addWidget(show_legend_,3,1,AlignTop);
			//legend
			label = new QLabel("Map intensity cutoff: ",background);
			grid->addWidget(label,4,0,AlignTop);
			intensity_cutoff_ = new QComboBox(false, background, "read-only combobox");
			intensity_cutoff_->insertItem("None");
			intensity_cutoff_->insertItem("Noise Estimator");
			grid->addWidget(intensity_cutoff_,4,1,AlignTop);
			
			grid->setRowStretch (5,2);
			
			tab->addTab(background,"General");
		
			//-------------DB Tab-------------
			background = new QWidget(tab);
			grid = new QGridLayout(background,3,3);
			grid->setMargin(6);
			grid->setSpacing(4);
		
			label = new QLabel("Host:",background);
			grid->addWidget(label,0,0,AlignLeft);
			db_host_ = new QLineEdit(background);
			grid->addWidget(db_host_,0,1,AlignLeft);
			db_host_->setMinimumWidth(fontMetrics().width('W') * 25);
		
			label = new QLabel("Port:",background);
			grid->addWidget(label,1,0,AlignLeft);
			db_port_ = new QLineEdit(background);
			grid->addWidget(db_port_,1,1,AlignLeft);
			db_port_->setMaximumWidth(fontMetrics().width('W') * 5);
		
			label = new QLabel("Database name:",background);
			grid->addWidget(label,2,0,AlignLeft);
			db_name_ = new QLineEdit(background);
			grid->addWidget(db_name_,2,1,AlignLeft);
		
			label = new QLabel("Login:",background);
			grid->addWidget(label,3,0,AlignLeft);
			db_login_ = new QLineEdit(background);
			grid->addWidget(db_login_,3,1,AlignLeft);
			
			grid->setRowStretch (4,2);
			
			tab->addTab(background,"DB");
		
			//-----------1D View Tab-----------
			background = new QWidget(tab);
			grid = new QGridLayout(background,2,2);
			grid->setMargin(6);
			grid->setSpacing(4);
		
			QGroupBox* box = new QGroupBox(2, Qt::Horizontal,"Colors",background);
			label = new QLabel("Peak color: ",box);
			peak_color_ = new ColorSelector(box);
			label = new QLabel("Icon color: ",box);
			icon_color_ = new ColorSelector(box);
			label = new QLabel("Highlighted peak color: ",box);
		
			high_color_ = new ColorSelector(box);
			label = new QLabel("Background color: ",box);
			back_color_1D_ = new ColorSelector(box);
			grid->addWidget(box,0,0);
		
			box = new QGroupBox(2,Qt::Horizontal,"Mapping",background);
			label = new QLabel("Map m/z to: ",box);
			axis_mapping_ = new QComboBox(false, box, "read-only combobox");
			axis_mapping_->insertItem("X-Axis");
			axis_mapping_->insertItem("Y-Axis");
			grid->addWidget(box,1,0);
			
			grid->setRowStretch (2,2);
			
			tab->addTab(background,"1D View");
		
			//-----------2D View Tab-----------
			background = new QWidget(tab);
		
			grid = new QGridLayout(background,2,2);
			grid->setMargin(6);
			grid->setSpacing(4);
			
			//dot mode
			box = new QGroupBox(2,Qt::Horizontal,"Dot coloring",background);
			QVButtonGroup* button_group;
			button_group = new QVButtonGroup("Mode:",box);
			button_group->setFrameStyle(QFrame::NoFrame);
			dot_mode_black_ = new QRadioButton("Black",button_group);
			dot_mode_gradient_ = new QRadioButton("Gradient",button_group);
			box->addSpace(0);
			dot_gradient_ = new MultiGradientSelector(box);
			grid->addMultiCellWidget(box,0,1,0,0);
		
			//surface mode
			box = new QGroupBox(2,Qt::Horizontal,"Surface coloring",background);
			surface_gradient_ = new MultiGradientSelector(box);
			grid->addWidget(box,0,1);
			
			//colors
			box = new QGroupBox(2,Qt::Horizontal,"Colors",background);
			label = new QLabel("Background color: ",box);
			back_color_2D_ = new ColorSelector(box);
			label = new QLabel("Interpolation steps: ",box);
			interpolation_steps_ = new QSpinBox(10,1000,1,box,"");
			interpolation_steps_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			grid->addWidget(box,1,1);

			//details
			box = new QGroupBox(2,Qt::Horizontal,"Surface/contour details",background);
			label = new QLabel("Squares per axis: ",box);
			marching_squares_steps_ = new QSpinBox(10,100,1,box,"");
			marching_squares_steps_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			label = new QLabel("Contour lines: ",box);
			contour_steps_ = new QSpinBox(3,30,1,box,"");
			contour_steps_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			grid->addWidget(box,2,0);
			
			//mapping
			box = new QGroupBox(2,Qt::Horizontal,"Mapping",background);
			label = new QLabel("Map m/z to: ",box);
			axis_mapping_2d_ = new QComboBox(false, box, "read-only combobox");
			axis_mapping_2d_->insertItem("X-Axis");
			axis_mapping_2d_->insertItem("Y-Axis");
			grid->addWidget(box,3,0);
			
			grid->setRowStretch (4,2);
			
			tab->addTab(background,"2D View");
	
			//-----------3D View Tab-----------
			background = new QWidget(tab);
		
			grid = new QGridLayout(background,4,2);
			grid->setMargin(6);
			grid->setSpacing(4);	
		
			box = new QGroupBox(2,Qt::Horizontal,"Dot coloring",background);
		  QVButtonGroup* coloring_group_3d = new QVButtonGroup("Color Mode:",box);
		  //coloring_group_3d->setFrameStyle(QFrame::NoFrame);
			box->addSpace(0);
	    
	    dot_mode_black_3d_ = new QRadioButton("Black",coloring_group_3d);
			dot_mode_gradient_3d_ = new QRadioButton("Gradient",coloring_group_3d);
			dot_gradient_3d_ = new MultiGradientSelector(coloring_group_3d);
			//box->addSpace(0);
			QHButtonGroup* interpolation_box = new QHButtonGroup("Interpolation steps",box);
			label = new QLabel("Interpolation steps: ",interpolation_box);
			dot_interpolation_steps_3d_ = new QSpinBox(10,1000,1,interpolation_box,"");
			dot_interpolation_steps_3d_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			box->addSpace(0);
			QVButtonGroup* shading_group_3d = new QVButtonGroup("Shade Mode:",box);
			//shading_group_3d->setFrameStyle(QFrame::NoFrame);
			shade_mode_flat_3d_ = new QRadioButton("Flat",shading_group_3d);
			shade_mode_smooth_3d_ = new QRadioButton("Smooth",shading_group_3d);	
			grid->addMultiCellWidget(box,0,3,0,0);
			
			box = new QGroupBox(2,Qt::Horizontal,"Line Width",background);
			label = new QLabel("Line Width: ",box);
			dot_line_width_ = new QSpinBox(1,10,1,box,"");
			dot_line_width_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			grid->addWidget(box,0,1);	
			
			box = new QGroupBox(2,Qt::Horizontal,"Colors",background);
			label = new QLabel("Background color: ",box);
			back_color_3d_ = new ColorSelector(box);
			label = new QLabel("Axes color: ",box);
			axes_color_3d_ = new ColorSelector(box);
			
			grid->addWidget(box,1,1);	


			box = new QGroupBox(2,Qt::Horizontal,"Data",background);
			QVButtonGroup* data_group = new QVButtonGroup("DataReduction",box);
			reduction_off_3d_ = new QRadioButton("Off",data_group);	
			reduction_on_max_3d_ = new QRadioButton("Max-Reduction",data_group);
			reduction_on_sum_3d_ = new QRadioButton("Sum-Reduction",data_group);
			box->addSpace(0);
			label = new QLabel("Number of peaks per Reductionstep:(Max-Red.) ",box);
			reduction_ratio_max_3d_ = new QSpinBox(10,100,1,box,"");

			label = new QLabel("m/z-Range per Reductionstep:(Sum-Red.)",box);
			reduction_ratio_sum_3d_ = new QSpinBox(10,100,1,box,"");
			grid->addWidget(box,2,1);	

			grid->addWidget(box,2,1);	
			grid->setRowStretch (2,2);
			
			tab->addTab(background,"3D View");
	
			load();
		}
	
		SpectrumMDIWindowPDP::~SpectrumMDIWindowPDP()
		{
		
		}
		
		void SpectrumMDIWindowPDP::load()
		{
			//general
			main_default_path_->setText(manager_->getPrefAsString("Preferences:DefaultPath").c_str());
			recent_files_->setValue(manager_->getPrefAsInt("Preferences:NumberOfRecentFiles"));
			default_map_view_->setCurrentText(manager_->getPrefAsString("Preferences:DefaultMapView").c_str());
			show_legend_->setCurrentText(manager_->getPrefAsString("Preferences:Legend").c_str());
			intensity_cutoff_->setCurrentText(manager_->getPrefAsString("Preferences:MapIntensityCutoff").c_str());
			
			
			//DB
			db_host_->setText(manager_->getPrefAsString("Preferences:DB:Host").c_str());
			db_port_->setText(manager_->getPrefAsString("Preferences:DB:Port").c_str());
			db_name_->setText(manager_->getPrefAsString("Preferences:DB:Name").c_str());
			db_login_->setText(manager_->getPrefAsString("Preferences:DB:Login").c_str());

			//1D
			peak_color_->setColor(QColor(manager_->getPrefAsString("Preferences:1D:PeakColor").c_str()));
			icon_color_->setColor(QColor(manager_->getPrefAsString("Preferences:1D:IconColor").c_str()));
			high_color_->setColor(QColor(manager_->getPrefAsString("Preferences:1D:HighColor").c_str()));
			back_color_1D_->setColor(QColor(manager_->getPrefAsString("Preferences:1D:BackgroundColor").c_str()));
			axis_mapping_->setCurrentText(manager_->getPrefAsString("Preferences:1D:Mapping:MappingOfMzTo").c_str());

			//2D
			if (UnsignedInt(manager_->getPref("Preferences:2D:Dot:Mode"))==Spectrum2DCanvas::DOT_GRADIENT)
			{
				dot_mode_gradient_->setChecked(true);
			}
			else if (UnsignedInt(manager_->getPref("Preferences:2D:Dot:Mode"))==Spectrum2DCanvas::DOT_BLACK)
			{
				dot_mode_black_->setChecked(true);
			}
			
			marching_squares_steps_->setValue(UnsignedInt(manager_->getPref("Preferences:2D:MarchingSquaresSteps")));
			contour_steps_->setValue(UnsignedInt(manager_->getPref("Preferences:2D:Contour:Lines")));
			dot_gradient_->gradient().fromString(manager_->getPrefAsString("Preferences:2D:Dot:Gradient"));
			interpolation_steps_->setValue(UnsignedInt(manager_->getPref("Preferences:2D:InterpolationSteps")));
			surface_gradient_->gradient().fromString(manager_->getPrefAsString("Preferences:2D:Surface:Gradient"));
			back_color_2D_->setColor(QColor(manager_->getPrefAsString("Preferences:2D:BackgroundColor").c_str()));
			axis_mapping_2d_->setCurrentText(manager_->getPrefAsString("Preferences:2D:Mapping:MappingOfMzTo").c_str());
	
			//3d

			if (UnsignedInt(manager_->getPref("Preferences:3D:Data:Mode"))==Spectrum3DCanvas::REDUCTION_MAX)
			{
					reduction_on_max_3d_->setChecked(true);
			}
			else if (UnsignedInt(manager_->getPref("Preferences:3D:Data:Mode"))==Spectrum3DCanvas::REDUCTION_OFF)
			{
				reduction_off_3d_->setChecked(true);
			}
			else if(UnsignedInt(manager_->getPref("Preferences:3D:Data:Mode"))==Spectrum3DCanvas::REDUCTION_SUM)
			{
				reduction_on_sum_3d_->setChecked(true);
			}
			
			reduction_ratio_max_3d_->setValue(UnsignedInt(manager_->getPref("Preferences:3D:Data:Reduction:Max")));
			reduction_ratio_sum_3d_->setValue(UnsignedInt(manager_->getPref("Preferences:3D:Data:Reduction:Sum")));
		
			if (UnsignedInt(manager_->getPref("Preferences:3D:Dot:Mode"))==Spectrum3DCanvas::DOT_GRADIENT)
			{
				dot_mode_gradient_3d_->setChecked(true);

				if (UnsignedInt(manager_->getPref("Preferences:3D:Shade:Mode"))==Spectrum3DCanvas::SHADE_FLAT)
				{
					shade_mode_flat_3d_->setChecked(true);
				}
				else if (UnsignedInt(manager_->getPref("Preferences:3D:Shade:Mode"))==Spectrum3DCanvas::SHADE_SMOOTH)
				{
					shade_mode_smooth_3d_->setChecked(true);
				}
			}
			else if (UnsignedInt(manager_->getPref("Preferences:3D:Dot:Mode"))==Spectrum3DCanvas::DOT_BLACK)
			{
				dot_mode_black_3d_->setChecked(true);
			}

			dot_interpolation_steps_3d_->setValue(UnsignedInt(manager_->getPref("Preferences:3D:Dot:InterpolationSteps")));
			back_color_3d_->setColor(QColor(manager_->getPrefAsString("Preferences:3D:BackgroundColor").c_str()));
			dot_gradient_3d_->gradient().fromString(manager_->getPrefAsString("Preferences:3D:Dot:Gradient"));
			axes_color_3d_->setColor(QColor(manager_->getPrefAsString("Preferences:3D:AxesColor").c_str()));
			dot_line_width_->setValue(UnsignedInt(manager_->getPref("Preferences:3D:Dot:LineWidth")));
		
		}
		
		void SpectrumMDIWindowPDP::save()
		{
			//main
			manager_->setPref("Preferences:DefaultPath", main_default_path_->text().ascii());
			manager_->setPref("Preferences:NumberOfRecentFiles", recent_files_->value());
			manager_->setPref("Preferences:DefaultMapView", default_map_view_->currentText().ascii());
			manager_->setPref("Preferences:Legend", show_legend_->currentText().ascii());
			manager_->setPref("Preferences:MapIntensityCutoff", intensity_cutoff_->currentText().ascii());

			//DB
			manager_->setPref("Preferences:DB:Host",db_host_->text().ascii());
			manager_->setPref("Preferences:DB:Port",db_port_->text().ascii());
			manager_->setPref("Preferences:DB:Name",db_name_->text().ascii());
			manager_->setPref("Preferences:DB:Login",db_login_->text().ascii());
			manager_->removePref("DBPassword");

			//1D
			manager_->setPref("Preferences:1D:PeakColor",peak_color_->getColor().name().ascii());
			manager_->setPref("Preferences:1D:IconColor",icon_color_->getColor().name().ascii());
			manager_->setPref("Preferences:1D:HighColor",high_color_->getColor().name().ascii());
			manager_->setPref("Preferences:1D:BackgroundColor",back_color_1D_->getColor().name().ascii());
			manager_->setPref("Preferences:1D:Mapping:MappingOfMzTo",axis_mapping_->currentText().ascii());
			
			//2D
			if (dot_mode_gradient_->isChecked())
			{
				manager_->setPref("Preferences:2D:Dot:Mode", Spectrum2DCanvas::DOT_GRADIENT);
			}
			else	if (dot_mode_black_->isChecked())
			{
				manager_->setPref("Preferences:2D:Dot:Mode", Spectrum2DCanvas::DOT_BLACK);
			}
			manager_->setPref("Preferences:2D:MarchingSquaresSteps",marching_squares_steps_->value());
			manager_->setPref("Preferences:2D:Contour:Lines",contour_steps_->value());
			manager_->setPref("Preferences:2D:InterpolationSteps",interpolation_steps_->value());
			manager_->setPref("Preferences:2D:BackgroundColor",back_color_2D_->getColor().name().ascii());
			manager_->setPref("Preferences:2D:Mapping:MappingOfMzTo",axis_mapping_2d_->currentText().ascii());
			manager_->setPref("Preferences:2D:Dot:Gradient",dot_gradient_->gradient().toString());
			manager_->setPref("Preferences:2D:Surface:Gradient",surface_gradient_->gradient().toString());
			
			//3d
			if (dot_mode_gradient_3d_->isChecked())
			{
				manager_->setPref("Preferences:3D:Dot:Mode", Spectrum3DCanvas::DOT_GRADIENT);
			}
			else	if (dot_mode_black_3d_->isChecked())
			{
				manager_->setPref("Preferences:3D:Dot:Mode", Spectrum3DCanvas::DOT_BLACK);	
			}
			manager_->setPref("Preferences:3D:Dot:Gradient",dot_gradient_3d_->gradient().toString());
			manager_->setPref("Preferences:3D:Dot:InterpolationSteps",dot_interpolation_steps_3d_->value());

			if (shade_mode_flat_3d_->isChecked())
			{
			manager_->setPref("Preferences:3D:Shade:Mode", Spectrum3DCanvas::SHADE_FLAT);
			}
			else	if (shade_mode_smooth_3d_->isChecked())
			{
				manager_->setPref("Preferences:3D:Shade:Mode", Spectrum3DCanvas::SHADE_SMOOTH);	
			}
		  manager_->setPref("Preferences:3D:BackgroundColor",back_color_3d_->getColor().name().ascii());
			manager_->setPref("Preferences:3D:AxesColor",axes_color_3d_->getColor().name().ascii());
			manager_->setPref("Preferences:3D:Dot:LineWidth",dot_line_width_->value());

			if(reduction_on_max_3d_->isChecked())
			{
				manager_->setPref("Preferences:3D:Data:Mode",	Spectrum3DCanvas::REDUCTION_MAX);
			}
			else if(reduction_off_3d_->isChecked())
			{
				manager_->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_OFF);
			}
			else if (reduction_on_sum_3d_->isChecked())
			{
				manager_->setPref("Preferences:3D:Data:Mode",	Spectrum3DCanvas::REDUCTION_SUM);
			}
			manager_->setPref("Preferences:3D:Data:Reduction:Max",reduction_ratio_max_3d_->value());
			manager_->setPref("Preferences:3D:Data:Reduction:Sum",reduction_ratio_sum_3d_->value());

		}
		
		void SpectrumMDIWindowPDP::browseDefaultPath_()
		{
			QString path = QFileDialog::getExistingDirectory(main_default_path_->text(), this, "get existing directory", "Choose a directory", TRUE );
			if (path!="")
			{
				main_default_path_->setText(path);
			}
		}

	} // namespace Internal

} //namespace


