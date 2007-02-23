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

#include <OpenMS/VISUAL/DIALOGS/TOPPViewBasePDP.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/ColorSelector.h>
// Qt
#include <QtGui/QLayout>
#include <QtGui/QLineEdit>
#include <QtGui/QRadioButton>
#include <QtGui/QTabWidget>
#include <QtGui/QLabel>
#include <QtGui/QGroupBox>
#include <QtGui/QCheckBox>
#include <QtGui/QFileDialog>
#include <QtGui/QComboBox>
#include <QtGui/QPushButton>
#include <QtGui/QSpinBox>
#include <QtGui/QGridLayout>

using namespace std;

namespace OpenMS
{

	namespace Internal
	{


		TOPPViewBasePDP::TOPPViewBasePDP( TOPPViewBase* manager, QWidget* parent)
			:PreferencesDialogPage(manager,parent)
		{
			help_ = "This is the preferences dialog of the main window!";
		
			//temporary box
			QGroupBox* box;
			
			//tab widget
			QTabWidget* tab_widget = new QTabWidget(this);
			(new QGridLayout(this))->addWidget(tab_widget);
			//tab	
			QWidget* tab; 
			//tab layout
			QGridLayout* grid;
			
			//-----------General Tab-----------
			tab = new QWidget(tab_widget);
			grid = new QGridLayout(tab);
			
			//default path
			main_default_path_ = new QLineEdit(tab);
			addWidget(grid,0,"Default path:",main_default_path_);
			main_default_path_->setMinimumWidth(fontMetrics().width('W') * 25);
			QPushButton* tmp = new QPushButton("Browse",tab);
			grid->addWidget(tmp,0,2);
			connect(tmp,SIGNAL(clicked()),this,SLOT(browseDefaultPath_()));
			//recent files
			recent_files_ = addSpinBox(tab,1,20,1);
			addWidget(grid,1,"Number of recent files:",recent_files_);
			//default map view
			default_map_view_ = new QComboBox( tab);
			default_map_view_->insertItem(0,"2D");
			default_map_view_->insertItem(1,"3D");
			addWidget(grid,2,"Default map visualization:",default_map_view_);
			//legend
			show_legend_ = new QComboBox( tab);
			show_legend_->insertItem(0,"Show");
			show_legend_->insertItem(1,"Hide");
			addWidget(grid,3,"Axis legend:",show_legend_);
			//legend
			intensity_cutoff_ = new QComboBox( tab);
			intensity_cutoff_->insertItem(0,"None");
			intensity_cutoff_->insertItem(1,"Noise Estimator");
			addWidget(grid,4,"Map intensity cutoff:",intensity_cutoff_);
			
			finish(grid);
			
			tab_widget->addTab(tab,"General");
		
			//-------------DB Tab-------------
			tab = new QWidget(tab_widget);
			grid = new QGridLayout(tab);
		
			db_host_ = new QLineEdit(tab);
			db_host_->setMaximumWidth(fontMetrics().width('W') * 15);
			addWidget(grid,0,"Host:",db_host_);
		
			db_port_ = new QLineEdit(tab);
			db_port_->setMaximumWidth(fontMetrics().width('W') * 5);
			addWidget(grid,1,"Port:",db_port_);
		
			db_name_ = new QLineEdit(tab);
			addWidget(grid,2,"Database name:",db_name_);
		
			db_login_ = new QLineEdit(tab);
			addWidget(grid,3,"Login:",db_login_);
			
			finish(grid);
			
			tab_widget->addTab(tab,"DB");
		
			//-----------1D View Tab-----------
			tab = new QWidget(tab_widget);
			grid = new QGridLayout(tab);
			
			//color box
			box = addBox(grid,0,0,"Colors");
			peak_color_ = new ColorSelector(box);
			addWidget(box->layout(),0,"Peak color:",peak_color_);
			icon_color_ = new ColorSelector(box);
			addWidget(box->layout(),1,"Icon color:",icon_color_);
			high_color_ = new ColorSelector(box);
			addWidget(box->layout(),2,"Highlighted peak color:",high_color_);
			back_color_1D_ = new ColorSelector(box);
			addWidget(box->layout(),3,"Background color:",back_color_1D_);						

			//mapping box
			box = addBox(grid,1,0,"Mapping");		
			axis_mapping_ = new QComboBox( box);
			axis_mapping_->insertItem(0,"X-Axis");
			axis_mapping_->insertItem(1,"Y-Axis");
			addWidget(box->layout(),0,"Map m/z to: ",axis_mapping_);
			
			finish(grid);
			
			tab_widget->addTab(tab,"1D View");
		
			//-----------2D View Tab-----------
			tab = new QWidget(tab_widget);
			grid = new QGridLayout(tab);

			//colors
			box = addBox(grid,0,0,"Colors");
			back_color_2D_ = new ColorSelector(box);
			addWidget(box->layout(),0,"Background color:",back_color_2D_);
			interpolation_steps_ = addSpinBox(box,10,1000,1);
			addWidget(box->layout(),1,"Interpolation steps:",interpolation_steps_);
			finish(box->layout());
			
			//mapping
			box = addBox(grid,1,0,"Mapping");
			axis_mapping_2d_ = new QComboBox( box);
			axis_mapping_2d_->insertItem(0,"X-Axis");
			axis_mapping_2d_->insertItem(1,"Y-Axis");
			addWidget(box->layout(),0,"Map m/z to:",axis_mapping_2d_);
			finish(box->layout());			
						
			//dot mode
			box = addBox(grid,0,1,"Dot Colors");
			QVBoxLayout* tmp2 = new QVBoxLayout();
			dot_mode_black_ = new QRadioButton("Black",tab);
			tmp2->addWidget(dot_mode_black_);
			dot_mode_gradient_ = new QRadioButton("Gradient",tab);
			tmp2->addWidget(dot_mode_gradient_);
			addLayout(box->layout(),0,"Mode:",tmp2);
	
			dot_gradient_ = new MultiGradientSelector(box);
			addWidget(box->layout(),1,"Gradient:",dot_gradient_);
			finish(box->layout());
					
			//surface mode
			box = addBox(grid,1,1,"Surface/contour settings");
			surface_gradient_ = new MultiGradientSelector(box);
			addWidget(box->layout(),0,"Gradient:",surface_gradient_);			

			marching_squares_steps_ = addSpinBox(box,10,100,1);
			addWidget(box->layout(),1,"Squares per axis:",marching_squares_steps_);
			
			contour_steps_ = addSpinBox(box,3,30,1);
			addWidget(box->layout(),2,"Contour lines:",contour_steps_);
			finish(box->layout());
			
			finish(grid);
			
			tab_widget->addTab(tab,"2D View");
	
			//-----------3D View Tab-----------
			tab = new QWidget(tab_widget);
			grid = new QGridLayout(tab);
			
			//peak color box
			box = addBox(grid,0,0,"Peak colors",1,2);

			tmp2 = new QVBoxLayout();
			dot_mode_black_3d_ = new QRadioButton("Black",tab);
			tmp2->addWidget(dot_mode_black_3d_);
			dot_mode_gradient_3d_ = new QRadioButton("Gradient",tab);
			tmp2->addWidget(dot_mode_gradient_3d_);
			addLayout(box->layout(),0,"Mode:",tmp2);

			dot_gradient_3d_ = new MultiGradientSelector(box);
			addWidget(box->layout(),1,"Gradient:",dot_gradient_3d_);

			dot_interpolation_steps_3d_ = addSpinBox(box,10,1000,1);
			addWidget(box->layout(),2,"Interpolation steps:",dot_interpolation_steps_3d_);
			finish(box->layout());
			
			tmp2 = new QVBoxLayout();
			QButtonGroup* tmp3 = new QButtonGroup(tab);
			shade_mode_flat_3d_ = new QRadioButton("Flat",tab);
			tmp2->addWidget(shade_mode_flat_3d_);
			tmp3->addButton(shade_mode_flat_3d_);
			shade_mode_smooth_3d_ = new QRadioButton("Smooth",tab);
			tmp2->addWidget(shade_mode_smooth_3d_);
			tmp3->addButton(shade_mode_smooth_3d_);
			addLayout(box->layout(),3,"Shade mode:",tmp2);
			finish(box->layout());
			
			//misc box
			box = addBox(grid,1,0,"Misc");
			
			dot_line_width_ = addSpinBox(box,1,10,1);
			addWidget(box->layout(),0,"Line width:",dot_line_width_);
			
			back_color_3d_ = new ColorSelector(box);
			addWidget(box->layout(),1,"Background color:",back_color_3d_);
			axes_color_3d_ = new ColorSelector(box);
			addWidget(box->layout(),2,"Axis color:",axes_color_3d_);
			finish(box->layout());

			//data reduction box
		 	box = addBox(grid,1,1,"Data reduction");
			data_reduction_3d_ = new QComboBox( box);
			data_reduction_3d_->insertItem(0,"Reduction OFF");
			data_reduction_3d_->insertItem(1,"MaxReduction");
			data_reduction_3d_->insertItem(2,"SumReduction");
			addWidget(box->layout(),0,"Mode:",data_reduction_3d_);

			reduction_diplay_peaks_3d_ = addSpinBox(box,5000,200000,5000);
			addWidget(box->layout(),1,"Displayed Peaks:",reduction_diplay_peaks_3d_);
			finish(box->layout());
			
			finish(grid);
						
			tab_widget->addTab(tab,"3D View");
	
			load();
		}
	
		TOPPViewBasePDP::~TOPPViewBasePDP()
		{
		
		}
		
		void TOPPViewBasePDP::load()
		{
			//general
			main_default_path_->setText(manager_->getPrefAsString("Preferences:DefaultPath").c_str());
			recent_files_->setValue(manager_->getPrefAsInt("Preferences:NumberOfRecentFiles"));
			default_map_view_->setCurrentIndex(default_map_view_->findText(manager_->getPrefAsString("Preferences:DefaultMapView").c_str()));
			show_legend_->setCurrentIndex(show_legend_->findText(manager_->getPrefAsString("Preferences:Legend").c_str()));
			intensity_cutoff_->setCurrentIndex(intensity_cutoff_->findText(manager_->getPrefAsString("Preferences:MapIntensityCutoff").c_str()));
			
			
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
			axis_mapping_->setCurrentIndex(axis_mapping_->findText(manager_->getPrefAsString("Preferences:1D:Mapping:MappingOfMzTo").c_str()));

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
			axis_mapping_2d_->setCurrentIndex(axis_mapping_2d_->findText(manager_->getPrefAsString("Preferences:2D:Mapping:MappingOfMzTo").c_str()));
	
			//3d
		
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
			
			data_reduction_3d_->setCurrentIndex(data_reduction_3d_->findText(manager_->getPrefAsString("Preferences:3D:Reduction:Mode").c_str()));	
			reduction_diplay_peaks_3d_->setValue(UnsignedInt(manager_->getPrefAsInt("Preferences:3D:DisplayedPeaks")));
		
			dot_interpolation_steps_3d_->setValue(UnsignedInt(manager_->getPref("Preferences:3D:Dot:InterpolationSteps")));
			back_color_3d_->setColor(QColor(manager_->getPrefAsString("Preferences:3D:BackgroundColor").c_str()));
			dot_gradient_3d_->gradient().fromString(manager_->getPrefAsString("Preferences:3D:Dot:Gradient"));
			axes_color_3d_->setColor(QColor(manager_->getPrefAsString("Preferences:3D:AxesColor").c_str()));
			dot_line_width_->setValue(UnsignedInt(manager_->getPref("Preferences:3D:Dot:LineWidth")));
		
		}
		
		void TOPPViewBasePDP::save()
		{
			//main
			manager_->setPref("Preferences:DefaultPath", main_default_path_->text().toAscii().data());
			manager_->setPref("Preferences:NumberOfRecentFiles", recent_files_->value());
			manager_->setPref("Preferences:DefaultMapView", default_map_view_->currentText().toAscii().data());
			manager_->setPref("Preferences:Legend", show_legend_->currentText().toAscii().data());
			manager_->setPref("Preferences:MapIntensityCutoff", intensity_cutoff_->currentText().toAscii().data());

			//DB
			manager_->setPref("Preferences:DB:Host",db_host_->text().toAscii().data());
			manager_->setPref("Preferences:DB:Port",db_port_->text().toAscii().data());
			manager_->setPref("Preferences:DB:Name",db_name_->text().toAscii().data());
			manager_->setPref("Preferences:DB:Login",db_login_->text().toAscii().data());
			manager_->removePref("DBPassword");

			//1D
			manager_->setPref("Preferences:1D:PeakColor",peak_color_->getColor().name().toAscii().data());
			manager_->setPref("Preferences:1D:IconColor",icon_color_->getColor().name().toAscii().data());
			manager_->setPref("Preferences:1D:HighColor",high_color_->getColor().name().toAscii().data());
			manager_->setPref("Preferences:1D:BackgroundColor",back_color_1D_->getColor().name().toAscii().data());
			manager_->setPref("Preferences:1D:Mapping:MappingOfMzTo",axis_mapping_->currentText().toAscii().data());
			
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
			manager_->setPref("Preferences:2D:BackgroundColor",back_color_2D_->getColor().name().toAscii().data());
			manager_->setPref("Preferences:2D:Mapping:MappingOfMzTo",axis_mapping_2d_->currentText().toAscii().data());
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
		  manager_->setPref("Preferences:3D:BackgroundColor",back_color_3d_->getColor().name().toAscii().data());
			manager_->setPref("Preferences:3D:AxesColor",axes_color_3d_->getColor().name().toAscii().data());
			manager_->setPref("Preferences:3D:Dot:LineWidth",dot_line_width_->value());
	
			if(data_reduction_3d_->currentText().toAscii().data()=="MaxReduction")
			{
				manager_->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_MAX);
			}
			else if(data_reduction_3d_->currentText().toAscii().data()=="Reduction OFF")
			{
				manager_->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_OFF);
			}
			else if(data_reduction_3d_->currentText().toAscii().data()=="SumReduction")
			{
				manager_->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_SUM);
			}
			
			manager_->setPref("Preferences:3D:Reduction:Mode", data_reduction_3d_->currentText().toAscii().data());
			manager_->setPref("Preferences:3D:DisplayedPeaks",	reduction_diplay_peaks_3d_->value());
		
		}
		
		void TOPPViewBasePDP::browseDefaultPath_()
		{
			QString path = QFileDialog::getExistingDirectory(this, "Choose a directory", main_default_path_->text());
			if (path!="")
			{
				main_default_path_->setText(path);
			}
		}

	} // namespace Internal

} //namespace


