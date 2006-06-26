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
// $Id: SpectrumMDIWindow.C,v 1.120 2006/06/09 22:00:08 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#ifdef ANDIMS_DEF
#include <OpenMS/FORMAT/ANDIFile.h>
#endif

#ifdef DB_DEF
#include <OpenMS/FORMAT/DBAdapter.h>
#include <OpenMS/VISUAL/DIALOGS/DBSpectrumSelectorDialog.h>
#endif

#include <OpenMS/VISUAL/DIALOGS/SaveImageDialog.h>
#include <OpenMS/VISUAL/SpectrumMDIWindow.h>
#include <OpenMS/VISUAL/DIALOGS/SpectrumMDIWindowPDP.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/DIALOGS/FeaFiDialog.h>
#include <OpenMS/VISUAL/Spectrum1DWindow.h>
#include <OpenMS/VISUAL/Spectrum2DWindow.h>
#include <OpenMS/VISUAL/Spectrum3DWindow.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>
#include <OpenMS/VISUAL/LayerManager.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/DPeakPickerCWT.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>

// Qt
#include <qaction.h>
#include <qmainwindow.h>
#include <qtoolbutton.h>
#include <qworkspace.h>
#include <qfiledialog.h>
#include <qapplication.h>
#include <qprinter.h>
#include <qpainter.h>
#include <qpaintdevicemetrics.h>
#include <qmessagebox.h>
#include <qlabel.h>
#include <qtooltip.h>
#include <qinputdialog.h>
#include <qstring.h>
#include <qmessagebox.h>
#include <qlineedit.h>
#include <qdir.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qstatusbar.h>
#include <qimage.h>


// STL
#include <iostream>
#include <sstream>

//xpm buttons
#include "ICONS/zoom.xpm"
#include "ICONS/translate.xpm"
#include "ICONS/noaction.xpm"
#include "ICONS/grid.xpm"
#include "ICONS/print.xpm"
#include "ICONS/colors.xpm"
#include "ICONS/contours.xpm"
#include "ICONS/lines.xpm"
#include "ICONS/peaks.xpm"
#include "ICONS/points.xpm"
#include "ICONS/intensity_scaled_dots.xpm"
#include "ICONS/reset_zoom.xpm"
#include "ICONS/tile_horizontal.xpm"
#include "ICONS/tile_vertical.xpm"
#include "ICONS/measure.xpm"
#include "ICONS/top_peaks.xpm"
#include "ICONS/3d_peaks.xpm"


using namespace std;

namespace OpenMS
{
	using namespace Internal;
	
	SpectrumMDIWindow* SpectrumMDIWindow::instance_ = 0;
	
	SpectrumMDIWindow* SpectrumMDIWindow::instance()
	{
		if (!instance_)
		{
			instance_ = new SpectrumMDIWindow();
		}
		return instance_;
	}
	
	SpectrumMDIWindow::SpectrumMDIWindow(QWidget* parent, const char* name, WFlags f):
		QMainWindow(parent,name,f),

		PreferencesManager(),
		recent_files_()
	{
		//prevents errors caused by too small width,height values
		setMinimumSize(200,200);
	
		// create dummy widget (to be able to have a layout), Tab bar and workspace
		QWidget* dummy = new QWidget(this);
		setCentralWidget(dummy);
		QVBoxLayout* box_layout = new QVBoxLayout(dummy);
		tab_bar_ = new EnhancedTabBar(dummy);
		QTab* tmp = new QTab("dummy");
		tab_bar_->addTab(tmp);
		tab_bar_->setMinimumSize(tab_bar_->sizeHint());
		tab_bar_->removeTab(tmp);
	
		//connect slots and sigals for selecting spectra
		connect(tab_bar_,SIGNAL(selected(int)),this,SLOT(focusSpectrumByAddress(int)));
		connect(tab_bar_,SIGNAL(doubleClicked(OpenMS::SignedInt)),this,SLOT(closeFileByTab(OpenMS::SignedInt)));
	
		box_layout->addWidget(tab_bar_);
		ws_=new QWorkspace(dummy);
		connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateToolbar(QWidget*)));
		connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateTabBar(QWidget*)));
		box_layout->addWidget(ws_);
	
		// File menu
		QPopupMenu * file = new QPopupMenu(this);
		menuBar()->insertItem("&File", file);
		file->insertItem("&Open",this,SLOT(openSpectrumDialog()));
		file->insertItem("&Close",this,SLOT(closeFile()));
		file->insertSeparator();
		recent_menu_ = new QPopupMenu(this);
		file->insertItem("Recent files",recent_menu_);
		file->insertSeparator();
		file->insertItem("&Preferences",this, SLOT(preferencesDialog()));
		file->insertItem("&Quit",qApp,SLOT(quit()));
	  //View menu
		QPopupMenu * view = new QPopupMenu(this);
		menuBar()->insertItem("&View", view);
		view->insertItem("&go to",this,SLOT(gotoDialog()), CTRL+Key_G);
	
		//Image menu
		QPopupMenu * image = new QPopupMenu(this);
		menuBar()->insertItem("&Image", image);
		image->insertItem("&Save to file",this,SLOT(saveImage()));
		image->insertItem("&Print",this,SLOT(print()));
	
		//Windows menu
		QPopupMenu * windows = new QPopupMenu(this);
		menuBar()->insertItem("&Windows", windows);
		windows->insertItem("&Cascade",this->ws_,SLOT(cascade()));
		windows->insertItem("&Tile automatic",this->ws_,SLOT(tile()));
		windows->insertItem(QIconSet(QPixmap(XPM_tile_h)),"&Tile vertical",this,SLOT(tileHorizontal()));
		windows->insertItem(QIconSet(QPixmap(XPM_tile_v)),"&Tile horizontal",this,SLOT(tileVertical()));
		//Tools menu
		tools_menu_ = new QPopupMenu(this);
		menuBar()->insertItem("&Tools", tools_menu_);
		tools_menu_->insertItem("&Show Selected Peaks of Active Spectrum as Peaklist", this, SLOT(showPeaklistActiveSpectrum()));
		tools_menu_->insertItem("&Pick Peaks", this, SLOT(pickActiveSpectrum()));
		tools_menu_->insertItem("&Find Features", this, SLOT(findFeaturesActiveSpectrum()));
	
		//create status bar
		message_label_ = new QLabel(statusBar());
		statusBar()->addWidget(message_label_,1);
	
		rt_label_ = new QLabel("RT: 12345678", statusBar());
		rt_label_->setMinimumSize(rt_label_->sizeHint());
		rt_label_->setText("");
		statusBar()->addWidget(rt_label_,0,true);
		mz_label_ = new QLabel("m/z: 12345678", statusBar());
		mz_label_->setMinimumSize(mz_label_->sizeHint());
		mz_label_->setText("");
		statusBar()->addWidget(mz_label_,0,true);
		int_label_ = new QLabel("Int: 123456789012", statusBar());
		int_label_->setMinimumSize(int_label_->sizeHint());
		int_label_->setText("");
		statusBar()->addWidget(int_label_,0,true);
	
		//create 1D + 2D toolbars and connect
		createToolBar_();
		connect(action_modes_,SIGNAL(selected(QAction*)),this,SLOT(setActionMode(QAction*)));
		connect(draw_modes_,SIGNAL(selected(QAction*)),this,SLOT(setDrawMode(QAction*)));
		connect(grid_button_,SIGNAL(toggled(bool)),this,SLOT(showGridLines(bool)));
		connect(grid_button_2d_,SIGNAL(toggled(bool)),this,SLOT(showGridLines(bool)));
		connect(action_modes_2d_,SIGNAL(selected(QAction*)),this,SLOT(setActionMode(QAction*)));
		connect(action_modes_3d_,SIGNAL(selected(QAction*)),this,SLOT(setActionMode(QAction*)));
	
		//layer bar (for managing several Spectra in a 1D-Window)
		layer_bar_ = new QToolBar(this,"layerbar");
		layer_bar_->setLabel( "Layer manager" );
		layer_bar_->setFrameStyle( QFrame::ToolBarPanel | QFrame::Plain );
	  layer_bar_->setLineWidth( 1 );
	  layer_bar_->setResizeEnabled(true);
		moveDockWindow(layer_bar_,DockRight);
		setDockEnabled(layer_bar_,DockTop,false);
		setDockEnabled(layer_bar_,DockBottom,false);
	
		//add Layer Manager
		layer_manager_ = new LayerManager(layer_bar_, "LayerManager");
		layer_bar_->hide();
		
		//register meta value names
		MetaInfo::registry().registerName("FeatureDrawMode", "Specify what to draw of the Feature: BoundingBox, ConvexHulls");

		//set preferences file name + load preferencs
		loadPreferences();
	}
	
	
	void SpectrumMDIWindow::addDBSpectrum(UnsignedInt db_id, bool as_new_window, bool maps_as_2d, bool maximize, OpenDialog::Mower /*use_mower*/)
	{
	#ifdef DB_DEF
	  QApplication::setOverrideCursor(Qt::WaitCursor);
		
		//DBConnection for all DB queries
		DBConnection con;
		con.connect(getPref("Preferences:DB:Name"), getPref("Preferences:DB:Login"),getPref("DBPassword"),getPref("Preferences:DB:Host"),getPrefAsInt("Preferences:DB:Port"));

		//DB adapter
		DBAdapter dba(con);	
		
		String db_id_string(db_id);
		con.executeQuery("SELECT count(id) from Spectrum where MSExperiment_id='"+db_id_string+"' and MS_Level='1'");
		
		//tab caption
		String caption = "DB ("+db_id_string+")";
		con.lastResult().first();
	
		SpectrumWindow* w;
		//open in new window
		if (as_new_window)
		{
			//create 1D View
			if (con.lastResult().value(0).toInt()==1)
			{
				// create 1D window
				w = new Spectrum1DWindow(ws_,"Spectrum1DWindow",WDestructiveClose);
				w->widget()->setMainPreferences(prefs_);
				Spectrum1DWindow* w1 = dynamic_cast<Spectrum1DWindow*>(w);
				
				//determine Spectrum id
				con.executeQuery("SELECT id from Spectrum where MSExperiment_id='"+db_id_string+"' and MS_Level='1'");
				con.lastResult().first();
				UID spectrum_id = con.lastResult().value(0).toInt();
				
				//improve the interface (template function) and hand it to canvas directly -> avoid copying the data
				Spectrum1DCanvas::ExperimentType& exp = w1->widget()->canvas()->addEmptyDataSet();
				exp.setName(caption);  // set layername
				MSSpectrum<>* spec = dba.loadSpectrum(spectrum_id);
				exp.push_back(*spec);
				delete(spec);
				w1->widget()->canvas()->finishAdding();
			}
			//create 2D/3D view
			else
			{
				//create 2D view
				if (maps_as_2d)
				{
					//create 2D window
					w = new Spectrum2DWindow(ws_,"Spectrum2DWindow",WDestructiveClose);
					w->widget()->setMainPreferences(prefs_);
					Spectrum2DWindow* w2 = dynamic_cast<Spectrum2DWindow*>(w);
					w2->widget()->canvas()->setDotMode(getPrefAsInt("Preferences:2D:Dot:Mode"));
					w2->widget()->canvas()->setDotGradient(getPrefAsString("Preferences:2D:Dot:Gradient"));
					w2->widget()->canvas()->setSurfaceGradient(getPref("Preferences:2D:Surface:Gradient"));
					
					//load spectrum
					MSExperiment<>* exp = dba.loadMSExperiment(db_id);
					
					//TODO remove: copy data
					Spectrum2DCanvas::ExperimentType& exp2 = w2->widget()->canvas()->addEmptyDataSet();
					exp2 = *exp;
					delete(exp);
					exp2.setName(caption);   // set layername
					w2->widget()->canvas()->finishAdding();
				}
				//create 3D view
				else
				{
					// create 3D window
					w = new Spectrum3DWindow(ws_, "Spectrum3DWindow", WDestructiveClose);
					w->widget()->setMainPreferences(prefs_);
					Spectrum3DWindow* w3 = dynamic_cast<Spectrum3DWindow*>(w);				
					//			w3->widget()->canvas()->setDotMode(getPrefAsInt("Preferences:3D:Dot:Mode"));
					w3->widget()->canvas()->setDotGradient(getPrefAsString("Preferences:3D:Dot:Gradient").c_str());
					connect(w3->widget()->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(update3DToolbar(QWidget*)));
		
					//load spectrum
					MSExperiment<>* exp = dba.loadMSExperiment(db_id);
					
					//TODO remove: copy data
					Spectrum3DCanvas::ExperimentType& exp2 = w3->widget()->canvas()->addEmptyDataSet();
					exp2 = *exp;
					delete(exp);
					exp2.setName(caption);   // set layername
					w3->widget()->canvas()->finishAdding();
				}
			}
		
			//do for all windows
			connectWindowSignals(w);
			w->setCaption(caption.c_str());
			addClient(w,caption);
			addTab(w,caption);
			
			//do for all (in active and in new window, 1D/2D/3D)
			if (w!=0)
			{
				if(maximize)
				{
					w->showMaximized();
				}
			}
		}
		//open in active window
		else
		{
			//create 1D View
			if (con.lastResult().value(0).toInt()==1)
			{
				Spectrum1DWindow* w1 = active1DWindow();
				//wrong active window type
				if (w1==0)
				{
					QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+db_id_string+") in a 2D window!<BR>Please open the file in new tab.").c_str());
				}
				else //open it
				{
					//determine Spectrum id
					con.executeQuery("SELECT id from Spectrum where MSExperiment_id='"+db_id_string+"' and MS_Level='1'");
					con.lastResult().first();
					UID spectrum_id = con.lastResult().value(0).toInt();
					
					//improve the interface (template function) and hand it to canvas directly -> avoid copying the data
					Spectrum1DCanvas::ExperimentType& exp = w1->widget()->canvas()->addEmptyDataSet();
					exp.setName(caption);  // set layername
					MSSpectrum<>* spec = dba.loadSpectrum(spectrum_id);
					exp.push_back(*spec);
					delete(spec);
					w1->widget()->canvas()->finishAdding();
				}
			}
			//create 2D/3D view
			else
			{
				Spectrum2DWindow* w2 = active2DWindow();
				Spectrum3DWindow* w3 = active3DWindow();
				//wrong active window type
				if (w2==0 && w3==0)
				{
					QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+db_id_string+") in a 2D/3D window!<BR>Please open the file in new tab.").c_str());
				}
				//create 2D view
				if (active2DWindow()!=0)
				{
					//load spectrum
					MSExperiment<>* exp = dba.loadMSExperiment(db_id);

					//Copy data
					Spectrum2DCanvas::ExperimentType& exp2 = w2->widget()->canvas()->addEmptyDataSet();
					exp2 = *exp;
					delete(exp);
					exp2.setName(caption);   // set layername
					w2->widget()->canvas()->finishAdding();
				}
				//create 3D view
				else
				{
					//load spectrum
					MSExperiment<>* exp = dba.loadMSExperiment(db_id);
	
					//TODO remove: copy data
					Spectrum3DCanvas::ExperimentType& exp2 = w3->widget()->canvas()->addEmptyDataSet();
					exp2 = *exp;
					delete(exp);
					exp2.setName(caption);   // set layername
					w3->widget()->canvas()->finishAdding();
				}
			}

			//do for all windows
			updateLayerbar();
		}
	
	  QApplication::restoreOverrideCursor();
	
	#endif
	}
	
	void SpectrumMDIWindow::preferencesDialog()
	{
		setActive(true);
		if (showPreferencesDialog())
		{
			savePreferences();
		}
	}
	
	void SpectrumMDIWindow::addSpectrum(const String& filename,bool as_new_window, bool maps_as_2d, bool maximize, OpenDialog::Mower use_mower)
	{
		//for comparison
		String filename_lower(filename);
		filename_lower.toLower();
		
	  setCursor(Qt::WaitCursor);
		if (filename != "")
		{
			QFileInfo file(filename.c_str());
			if (file.exists())
			{
				//extract the filename without path
				String caption = file.fileName().ascii();
				
				//update recent files list
				addRecentFile_(filename);
	
				//windowpointer
				SpectrumWindow* w=0;
	
				// 1D Files
				if (filename_lower.hasSuffix(".dta"))
				{
					//wrong active window type
					if (!as_new_window && active1DWindow()==0)
					{
						QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+filename+") in a 2D window!<BR>Please open the file in new tab.").c_str());
					}
					else //open it
					{
						w = addSpectrum1D_(filename,caption,as_new_window,use_mower);
					}
				}
				// 2D Files
				else
				{
					//open in 2D view
					if ((as_new_window && maps_as_2d) || (!as_new_window && active2DWindow()!=0))
					{
						w = addSpectrum2D_(filename,caption,as_new_window,use_mower);
					}
					//open in 3D view
					else if ((as_new_window && !maps_as_2d) || (!as_new_window && active3DWindow()!=0))
					{
						w = addSpectrum3D_(filename,caption,as_new_window,use_mower);
					}
					else
					//wrong active window type
					{
						QMessageBox::warning(this,"Wrong file type",("You cannot open 2D data ("+filename+") in a 1D window!<BR>Please open the file in new tab.").c_str());
					}
				}
				//do for all (in active and in new window, 1D/2D/3D)
				if (w!=0)
				{
					if(maximize)
					{
						w->showMaximized();
					}
				}
			}
			else //file does not exist
			{
				QMessageBox::warning(this,"Command line file",("The file "+filename+" does not exist!").c_str());
			}
		}
	  setCursor(Qt::ArrowCursor);
	}
	
	SpectrumWindow* SpectrumMDIWindow::addSpectrum1D_(const String& filename, const String& caption, bool as_new_window, OpenDialog::Mower /*use_mower*/)
	{
		//for comparison
		String filename_lower(filename);
		filename_lower.toLower();
		
		Spectrum1DWindow* w1;
		
		//open in active window
		if (as_new_window == false)
		{
			w1 = active1DWindow();
	
			//DTA
			if (filename_lower.hasSuffix(".dta"))
			{
				try
				{
					Spectrum1DCanvas::ExperimentType& exp = w1->widget()->canvas()->addEmptyDataSet();
					exp.setName(caption);  // set layername
					exp.resize(1);
					DTAFile inF;
					inF.load(filename,exp[0]);
					w1->widget()->canvas()->finishAdding();
					updateLayerbar();
				}
				catch(Exception::Base& e)
				{
					//QMessageBox::warning ( this , "Error while reading DTA file", e.what());
					cout << "Error while reading DTA file: " <<e.what()<<endl;
					return 0;
				}
			}
		}
		//open in new window
		else
		{
			w1 = new Spectrum1DWindow(ws_,"Spectrum1DWindow",WDestructiveClose);
			w1->widget()->setMainPreferences(prefs_);
			//try to read the data from file
			try
			{			
				//DTA
				if (filename_lower.hasSuffix(".dta"))
				{
					Spectrum1DCanvas::ExperimentType& exp = w1->widget()->canvas()->addEmptyDataSet();
					exp.setName(caption);  // set layername
					exp.resize(1);
					DTAFile inF;
					inF.load(filename,exp[0]);
					w1->widget()->canvas()->finishAdding();
					updateLayerbar();
				}
			}
			catch(Exception::Base& e)
			{
				//QMessageBox::warning ( this , "Error while reading file", e.what());
				cout << "Error while reading 1D file: " << e.what()<<endl;
				return 0;
			}
			
			connect(w1->widget()->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(update1DToolbar(QWidget*)));
						
			connectWindowSignals(w1);
			w1->setCaption(filename.c_str());
			addClient(w1,filename);
			addTab(w1,caption);
		}
	
		return w1;
	}
	
	SpectrumWindow* SpectrumMDIWindow::addSpectrum2D_(const String& filename, const String& caption, bool as_new_window, OpenDialog::Mower /*use_mower*/)
	{
		//for comparison
		String filename_lower(filename);
		filename_lower.toLower();
		
		Spectrum2DWindow* w2;
	
		//open in active window
		if (as_new_window == false)
		{
			w2 = active2DWindow();
	
			try
			{
				// FeatureMapFile
				if ((filename_lower.hasSuffix(".feat")))
				{
					DFeatureMap<2> map;
					DFeatureMapFile().load(filename,map);
					Spectrum2DCanvas::ExperimentType exp;
					map.sortByPosition();
					exp.set2DData(map);
					setFeatureMap_(w2->widget()->canvas(), exp, caption);
				}
				else
				{
					//do for each kind of 2D window
					Spectrum2DCanvas::ExperimentType& exp = w2->widget()->canvas()->addEmptyDataSet();
	
					//DTA2D
					if (filename_lower.hasSuffix(".dta2d"))
					{
						DTA2DFile().load(filename,exp);
					}
					//NETcdf
					else if ((filename_lower.hasSuffix(".cdf")))
					{
	#ifdef ANDIMS_DEF
						ANDIFile().load(filename,exp);
	#endif
					}
					//mzXML
					else if ((filename_lower.hasSuffix(".mzxml")))
					{
						MzXMLFile().load(filename,exp);
					}
					//mzData
					else if ((filename_lower.hasSuffix(".mzdata")))
					{
						MzDataFile().load(filename,exp);
					}
	
					//do for each kind of 2D window
					exp.setName(caption);  // set layername
					w2->widget()->canvas()->finishAdding();
				}
				updateLayerbar();
			}
			catch(Exception::Base& e)
			{
				//QMessageBox::warning ( this , "Error while reading DTA file", e.what());
				cout << "Error while reading 2D file: " << e.what()<<endl;
				return 0;
			}
		}
		//open in new window
		else
		{
			w2 = new Spectrum2DWindow(ws_, "Spectrum2DWindow", WDestructiveClose);
			
			w2->widget()->setMainPreferences(prefs_);					
			//try to read the data from file
			try
			{
				if ((filename_lower.hasSuffix(".feat")))
				{
					DFeatureMap<2> map;
					DFeatureMapFile().load(filename,map);
					Spectrum2DCanvas::ExperimentType exp;
					map.sortByPosition();
					exp.set2DData(map);
					setFeatureMap_(w2->widget()->canvas(), exp, caption);
				}
				else
				{
					Spectrum2DCanvas::ExperimentType& exp = w2->widget()->canvas()->addEmptyDataSet();
	
					//DTA2D
					if (filename_lower.hasSuffix(".dta2d"))
					{
						DTA2DFile().load(filename,exp);
					}
					//NETcdf
					else if ((filename_lower.hasSuffix(".cdf")))
					{
	#ifdef ANDIMS_DEF
						ANDIFile().load(filename,exp);
	#endif
					}
					//mzXML
					else if ((filename_lower.hasSuffix(".mzxml")))
					{
						MzXMLFile().load(filename,exp);
					}
					//mzData
					else if ((filename_lower.hasSuffix(".mzdata")))
					{
						MzDataFile().load(filename,exp);
					}
					
					// Do for each kind of file
					exp.setName(caption);  // set layername
					w2->widget()->canvas()->finishAdding();
				}
			}
			catch(Exception::Base& e)
			{
				//QMessageBox::warning ( this , "Error while reading file", e.what());
				cout << "Error while reading 2D file: " <<e.what()<<endl;
				return 0;
			}
						
			w2->widget()->canvas()->setDotMode(getPrefAsInt("Preferences:2D:Dot:Mode"));
			w2->widget()->canvas()->setDotGradient(getPrefAsString("Preferences:2D:Dot:Gradient"));
			w2->widget()->canvas()->setSurfaceGradient(getPrefAsString("Preferences:2D:Surface:Gradient"));
			connect(w2->widget()->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(update2DToolbar(QWidget*)));
				
			connectWindowSignals(w2);
			w2->setCaption(filename.c_str());
			addClient(w2,filename);			
			addTab(w2,caption);
		}	
		
		return w2;
	}
	
	SpectrumWindow* SpectrumMDIWindow::addSpectrum3D_(const String& filename, const String& caption, bool as_new_window, OpenDialog::Mower /*use_mower*/)
	{
		//for comparison
		String filename_lower(filename);
		filename_lower.toLower();
		
		Spectrum3DWindow* w3;
		
		//open in active window
		if (as_new_window == false)
		{
			w3 = active3DWindow();
			try
			{
				//do for each kind of 3D window
				Spectrum3DCanvas::ExperimentType& exp = w3->widget()->canvas()->addEmptyDataSet();
		
				//DTA2D
				if (filename_lower.hasSuffix(".dta2d"))
				{
					DTA2DFile().load(filename,exp);
				}
				//NETcdf
				else if ((filename_lower.hasSuffix(".cdf")))
				{
#ifdef ANDIMS_DEF
					ANDIFile().load(filename,exp);
#endif
				}
				//mzXML
				else if ((filename_lower.hasSuffix(".mzxml")))
				{
					MzXMLFile().load(filename,exp);
				}
				//mzData
				else if ((filename_lower.hasSuffix(".mzdata")))
				{
					MzDataFile().load(filename,exp);
				}
				
				//do for each kind of 3D window
				exp.setName(caption);  // set layername
				w3->widget()->canvas()->finishAdding();
				updateLayerbar();
			}
			catch(Exception::Base& e)
			{
				//QMessageBox::warning ( this , "Error while reading DTA file", e.what());
				cout << "Error while reading 2D file: " << e.what()<<endl;
				return 0;
			}
		}
		//open in new window
		else
		{
			w3 = new Spectrum3DWindow(ws_, "Spectrum3DWindow", WDestructiveClose);
			
			w3->widget()->setMainPreferences(prefs_);					
			//try to read the data from file
			try
			{
				Spectrum3DCanvas::ExperimentType& exp = w3->widget()->canvas()->addEmptyDataSet();
	
				//DTA2D
				if (filename_lower.hasSuffix(".dta2d"))
				{
					DTA2DFile().load(filename,exp);
				}
				//NETcdf
				else if ((filename_lower.hasSuffix(".cdf")))
				{
	#ifdef ANDIMS_DEF
					ANDIFile().load(filename,exp);
	#endif
				}
				//mzXML
				else if ((filename_lower.hasSuffix(".mzxml")))
				{
					MzXMLFile().load(filename,exp);
				}
				//mzData
				else if ((filename_lower.hasSuffix(".mzdata")))
				{
					MzDataFile().load(filename,exp);
				}
				
				exp.setName(caption);  // set layername
				w3->widget()->canvas()->finishAdding();
			}
			catch(Exception::Base& e)
			{
				//QMessageBox::warning ( this , "Error while reading file", e.what());
				cout << "Error while reading 2D file: " <<e.what()<<endl;
				return 0;
			}
			w3->widget()->canvas()->setDotGradient(getPrefAsString("Preferences:3D:Dot:Gradient"));
		
			connect(w3->widget()->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(update3DToolbar(QWidget*)));
			
			connectWindowSignals(w3);
			
			w3->setCaption(filename.c_str());
			
			addClient(w3,filename);
			addTab(w3,caption);
		}		
	
		return w3;
	}

	
	void SpectrumMDIWindow::addRecentFile_(const String& filename)
	{
		//get number of files in list from settings
		UnsignedInt number_of_recent_files_ = UnsignedInt(prefs_.getValue("Preferences:NumberOfRecentFiles"));
		String tmp = filename;
	
		//add prefix to relative paths
		if (!tmp.hasPrefix("/"))
		{
			tmp = QDir::currentDirPath().ascii()+string("/")+ tmp;
		}
	
		//check if the file is already in the vector. if so remove it
		recent_files_.erase(remove(recent_files_.begin(),recent_files_.end(),tmp),recent_files_.end());
	
		//add the file to the front
		recent_files_.insert(recent_files_.begin(),tmp);
	
		//shrink the recent_files to 10
		if (recent_files_.size()>number_of_recent_files_)
		{
			recent_files_.resize(number_of_recent_files_);
		}
	
		//update recent file menu
		updateRecentMenu_();
	}
	
	void SpectrumMDIWindow::updateRecentMenu_()
	{
		//update recent file menu
		recent_menu_->clear();
		for (UnsignedInt i=0; i<recent_files_.size(); ++i)
		{
			recent_menu_->insertItem(recent_files_[i].c_str(),this,SLOT(openRecentFile(int)),0,i);
		}
	}
	
	
	void SpectrumMDIWindow::addTab(SpectrumWindow* w, const String& tabCaption)
	{
		//add to tab bar (map the address of the widget and the pointer)
		QTab* tab = new QTab(tabCaption.c_str());
		tab_bar_->addTab(tab);
		tab->setIdentifier(PointerSizeInt(&(*w)));
		id_map_[PointerSizeInt(&(*w))]=w;
		//connect slots and sigals for removing the spectrum from the bar, when it is closed
		connect(w,SIGNAL(destroyed(QObject*)),this,SLOT(removeWidgetFromBar(QObject*)));
		tab_bar_->setCurrentTab(tab);
		
	}
	
	void SpectrumMDIWindow::maximizeActiveSpectrum()
	{
		if (ws_->activeWindow())
		{
			ws_->activeWindow()->showMaximized();
		}
	}
	
	void SpectrumMDIWindow::focusSpectrumByAddress(int address)
	{
	
		QWidget* m_w = id_map_[address];
		if (m_w)
		{
			m_w->setFocus();
		}
	}
	
	void SpectrumMDIWindow::saveImage()
	{
		//check if there is a active window
		SpectrumWindow* window = activeWindow();
		if (window!=0)
		{
			SaveImageDialog* dialog = new SaveImageDialog(this);
			dialog->setSize(1024,768);
			if (dialog->exec())
			{
				QString format=dialog->getFormat();
				QString file_name = QFileDialog::getSaveFileName( window->caption().section('.',0,0),
																													format+" Images(*."+format+" *."+format.lower()+")",
																													this, "save image dialog","Select file to save" );
				if (!file_name.isEmpty())
				{
					//append missing file extension
					if (!file_name.upper().endsWith(format))
					{
						file_name.append("."+format.lower());
					}
	
					QImage image = window->widget()->getImage(dialog->getXSize(),dialog->getYSize());
					image.save(file_name,format,100);
				}
			}
		}
	}
	
	void SpectrumMDIWindow::print()
	{
		#ifndef QT_NO_PRINTER
		SpectrumWindow* window = activeWindow();
		if (window!=0)
		{
			QPrinter* printer = new QPrinter(QPrinter::HighResolution);
			printer->setResolution(300);
			if (printer->setup(this))
			{
		 		QPainter p;
				if (!p.begin(printer)) return;
				QPaintDeviceMetrics metrics(p.device());
				unsigned int dpix = metrics.logicalDpiX();
				unsigned int dpiy = metrics.logicalDpiY();
				QRect body(dpix,dpiy,metrics.width()-2*dpix,metrics.height()-2*dpiy);			// one inch margin
				QImage image = window->widget()->getImage(body.width(),body.height());
				p.drawImage(body,image);
				QString titel = QString("%1\n%2").arg(window->caption().section('/',-1)).arg(QDate::currentDate().toString());
				//	p.drawText(dpix,0,body.width(),dpiy, Qt::AlignCenter, titel);
				p.drawText(dpix,body.height()+dpiy, body.width(), dpiy, Qt::AlignCenter, titel);
			}
			delete(printer);
		}
		#endif
	}
	
	void SpectrumMDIWindow::closeFile()
	{
		//check if there is a active window
		if (ws_->activeWindow())
		{
			ws_->activeWindow()->close();
		}
	}
	
	void SpectrumMDIWindow::closeFileByTab(SignedInt tab_identifier_)
	{
		SpectrumWindow* window = id_map_[tab_identifier_];
		window->close();
	}
	
	SpectrumMDIWindow::~SpectrumMDIWindow()
	{
		//do not add anything here. This class is a singleton
		//use closeEvent(...) to do the cleanup
	}
	
	void SpectrumMDIWindow::createToolBar_()
	{
		tool_bar_ = new QToolBar(this,"toolbar");
	
		action_modes_ = new QActionGroup(tool_bar_);
		action_modes_->setExclusive(TRUE);
	
		set_pick_action_ = new QAction( QString("Select"), QPixmap(XPM_noAction), NULL, CTRL + Key_Q, action_modes_,"SpectrumCanvas::AM_SELECT",TRUE);
		set_pick_action_->setOn(true);
		set_pick_action_->addTo(tool_bar_);
	
		set_zoom_action_ = new QAction( QString("Zoom"), QPixmap(XPM_zoom), NULL, CTRL + Key_W, action_modes_,"SpectrumCanvas::AM_ZOOM",TRUE);
		set_zoom_action_->addTo(tool_bar_);
	
		set_translate_action_ = new QAction( QString("Translate"), QPixmap(XPM_translate), NULL, CTRL + Key_R, action_modes_,"SpectrumCanvas::AM_TRANSLATE",TRUE);
		set_translate_action_->addTo(tool_bar_);
	
		tool_bar_->addSeparator();
	
		draw_modes_ = new QActionGroup(tool_bar_);
		draw_modes_->setExclusive(TRUE);
	
		set_peak_mode_ = new QAction( QString("Show peaks"), QPixmap(XPM_peaks), NULL, CTRL + Key_I, draw_modes_,"DM_PEAKS",TRUE);
		set_peak_mode_->addTo(tool_bar_);
	
		set_connected_lines_mode_ = new QAction( QString("Show connected lines"), QPixmap(XPM_lines), NULL, CTRL + Key_O, draw_modes_,"DM_CONNECTEDLINES",TRUE);
		set_connected_lines_mode_->addTo(tool_bar_);
	
		tool_bar_->addSeparator();
	
		grid_button_ = new QToolButton(QIconSet(QPixmap(XPM_grid)),"Show grid","Show grid",NULL,NULL,tool_bar_,"gridButton");
		grid_button_->setToggleButton(true);
		grid_button_->setOn(true);
	
		///	reset Zoom button
		tool_bar_->addSeparator();
		reset_zoom_button_ = new QToolButton(QIconSet(QPixmap(XPM_reset_zoom)),"Reset Zoom", "Reset Zoom", NULL, NULL, tool_bar_, "resetZoomButton");
		connect(reset_zoom_button_,SIGNAL(clicked()),this,SLOT(resetZoom()));
	
		tool_bar_->addSeparator();
	
		print_button_ = new QToolButton(QIconSet(QPixmap(XPM_print)),"Print","print",NULL,NULL,tool_bar_,"printButton");
		connect(print_button_,SIGNAL(clicked()),this,SLOT(print()));
	
		tool_bar_->addSeparator();
		link_box_ = new QComboBox(tool_bar_);
		QToolTip::add(link_box_,"Use this combobox to link two spectra.\nLinked spectra zoom in/out together");
		connect(link_box_,SIGNAL(activated(const QString&)),this,SLOT(linkActiveTo(const QString&)));
	
		tool_bar_->resize(tool_bar_->sizeHint());
		tool_bar_->show();
	
		// 2DWidget
		tool_bar_2d_ = new QToolBar(this, "toolbar2d");
	
		action_modes_2d_ = new QActionGroup(tool_bar_2d_);
		action_modes_2d_->setExclusive(TRUE);
	
		set_pick_action_2d_ = new QAction( QString("Select"), QPixmap(XPM_noAction), NULL, CTRL + Key_Q, action_modes_2d_,"SpectrumCanvas::AM_SELECT",TRUE);
		set_pick_action_2d_->setOn(true);
		set_pick_action_2d_->addTo(tool_bar_2d_);
	
		set_zoom_action_2d_ = new QAction( QString("Zoom"), QPixmap(XPM_zoom), NULL, CTRL + Key_W, action_modes_2d_,"SpectrumCanvas::AM_ZOOM",TRUE);
		set_zoom_action_2d_->addTo(tool_bar_2d_);
	
		set_translate_action_2d_ = new QAction( QString("Translate"), QPixmap(XPM_translate), NULL, CTRL + Key_R, action_modes_2d_,"SpectrumCanvas::AM_TRANSLATE",TRUE);
		set_translate_action_2d_->addTo(tool_bar_2d_);
	
		set_measure_action_2d_ = new QAction( QString("Measure"), QPixmap(XPM_measure), NULL, CTRL + Key_M, action_modes_2d_,"SpectrumCanvas::AM_MEASURE",TRUE);
		set_measure_action_2d_->addTo(tool_bar_2d_);
	
		tool_bar_2d_->addSeparator();
	
		show_points_button_2d_ = new QToolButton(QIconSet(QPixmap(XPM_points)), "Show dots", "Show dots", 0, 0, tool_bar_2d_, "showPoints");
		show_points_button_2d_->setToggleButton(true);
		show_points_button_2d_ ->setOn(false);
		connect(show_points_button_2d_, SIGNAL(toggled(bool)), this, SLOT(showPoints(bool)));
	
		intensity_scaled_dots_button_2d_ = new QToolButton(QIconSet(QPixmap(XPM_intensity_scaled_dots)), "Show intensisty scaled dots", "Show intensisty scaled dots", 0, 0, tool_bar_2d_, "setIntensityScaledDots");
		intensity_scaled_dots_button_2d_->setToggleButton(true);
		intensity_scaled_dots_button_2d_ ->setOn(false);
		connect(intensity_scaled_dots_button_2d_, SIGNAL(toggled(bool)), this, SLOT(setIntensityScaledDots(bool)));
	
		show_colors_button_2d_ = new QToolButton(QIconSet(QPixmap(XPM_colors)), "Show colored surface", "Show colored surface", 0, 0, tool_bar_2d_, "showColors");
		show_colors_button_2d_->setToggleButton(true);
		show_colors_button_2d_->setOn(false);
		connect(show_colors_button_2d_, SIGNAL(toggled(bool)), this, SLOT(showColors(bool)));
	
		show_contours_button_2d_ = new QToolButton(QIconSet(QPixmap(XPM_contours)), "Show contour lines", "Show contour lines", 0, 0, tool_bar_2d_, "showContours");
		show_contours_button_2d_->setToggleButton(true);
		show_contours_button_2d_->setOn(false);
		connect(show_contours_button_2d_, SIGNAL(toggled(bool)), this, SLOT(showContours(bool)));

		tool_bar_2d_->addSeparator();
	
		grid_button_2d_ = new QToolButton(QIconSet(QPixmap(XPM_grid)),"Show grid","Show grid",NULL,NULL,tool_bar_2d_,"gridButton");
		grid_button_2d_->setToggleButton(true);
		grid_button_2d_->setOn(true);
	
		tool_bar_2d_->addSeparator();
	
		reset_zoom_button_2d_ = new QToolButton(QIconSet(QPixmap(XPM_reset_zoom)), "Reset Zoom", "Reset Zoom", this, SLOT(resetZoom()), tool_bar_2d_, "resetZoom");
	

		//3Dwidget
		tool_bar_3d_ = new QToolBar(this, "toolbar3d");
		action_modes_3d_ = new QActionGroup(tool_bar_3d_);
		action_modes_3d_->setExclusive(TRUE);
		
		set_pick_action_3d_ = new QAction( QString("Select"), QPixmap(XPM_noAction), NULL, CTRL + Key_Q, action_modes_3d_,"SpectrumCanvas::AM_SELECT",TRUE);
		set_pick_action_3d_->setOn(true);
		set_pick_action_3d_->addTo(tool_bar_3d_);
		set_zoom_action_3d_ = new QAction( QString("Zoom"), QPixmap(XPM_zoom), NULL, CTRL + Key_W, action_modes_3d_,"SpectrumCanvas::AM_ZOOM",TRUE);
		set_zoom_action_3d_->addTo(tool_bar_3d_);
		tool_bar_3d_->addSeparator();
		
		show_back_view_3d_ =new QToolButton(QIconSet(QPixmap(XPM_3d_peaks)), "Side view", "Side view", 0, 0, tool_bar_3d_, "backView");
		show_back_view_3d_->setToggleButton(true);
		show_back_view_3d_->setOn(false);
		connect(show_back_view_3d_, SIGNAL(toggled(bool)), this, SLOT(setBackView3D(bool)));
			
		show_top_view_3d_=new QToolButton(QIconSet(QPixmap(XPM_top_peaks)), "Top view", "Top view", 0, 0, tool_bar_3d_, "showTopView");
		show_top_view_3d_->setToggleButton(true);
		show_top_view_3d_->setOn(false);		
		connect(show_top_view_3d_, SIGNAL(toggled(bool)), this, SLOT(setTopView3D(bool)));
		
		intensity_scaled_dots_button_3d_ = new QToolButton(QIconSet(QPixmap(XPM_intensity_scaled_dots)), "setIntensityScaledDots", "setIntensityScaledDots", 0, 0, tool_bar_3d_, "setIntensityScaledDots");
		intensity_scaled_dots_button_3d_->setToggleButton(true);
		intensity_scaled_dots_button_3d_ ->setOn(false);
		connect(intensity_scaled_dots_button_3d_, SIGNAL(toggled(bool)), this, SLOT(setIntensityScaledDots3D(bool)));
	
		show_reset_view_3d_ = new QToolButton(QIconSet(QPixmap(XPM_reset_zoom)), "Reset zoom", "Reset zoom", 0, 0, tool_bar_3d_, "resetZoom");
		connect(show_reset_view_3d_, SIGNAL(toggled(bool)), this, SLOT(resetZoom()));
		
}
	
	void SpectrumMDIWindow::linkActiveTo(const QString& path)
	{
		unlinkActive_();
		if (path!="<unlinked>")
		{
			QFileInfo file;
			//link spectras
			QWidgetList windows = ws_->windowList();
			for ( int i = 0; i < int(windows.count()); ++i )
			{
				QWidget *window = windows.at(i);
				file = QFileInfo(window->caption());
				if (file.fileName()+"  ("+window->caption()+")"==path)
				{
					//connect the slots
					connect(activeWindow()->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),dynamic_cast<SpectrumWindow*>(window)->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
					connect(dynamic_cast<SpectrumWindow*>(window)->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),activeWindow()->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
					//add links to the map
					link_map_[PointerSizeInt(&(*ws_->activeWindow()))]=PointerSizeInt(&(*window));
					link_map_[PointerSizeInt(&(*window))]=PointerSizeInt(&(*ws_->activeWindow()));
				}
			}
		}
	
	}
	
	void SpectrumMDIWindow::unlinkActive_()
	{
		int active_address = PointerSizeInt(&(*ws_->activeWindow()));
		int active_linked_to_address = link_map_[active_address];
		if (active_linked_to_address != 0)
		{		
			//remove signals		
			disconnect(id_map_[active_linked_to_address]->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)), activeWindow()->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
			disconnect(activeWindow()->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)), id_map_[active_linked_to_address]->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
			//remove from the map
			link_map_.erase(active_address);
			link_map_.erase(active_linked_to_address);
		}
	}
	
	void SpectrumMDIWindow::showStatusMessage(string msg,UnsignedInt time=0)
	{
		if (time==0)
		{
				message_label_->setText(msg.c_str());
				statusBar()->repaint();
		}
		else
		{
				statusBar()->message(msg.c_str(), time);
		}
	}
	
	void SpectrumMDIWindow::showCursorStatus(double mz, double intensity, double rt)
	{
		message_label_->setText("");
		if (mz==-1)
		{
			mz_label_->setText("m/z: ");
		}
		else
		{
			mz_label_->setText(("m/z: "+String(mz,8).fillLeft(' ',8)).c_str());
		}
		if (rt==-1)
		{
			rt_label_->setText("RT: ");
		}
		else
		{
			rt_label_->setText(("RT: "+String(rt,8).fillLeft(' ',8)).c_str());
		}
		if (intensity==-1)
		{
			int_label_->setText("Int: ");
		}
		else
		{
			int_label_->setText(("Int: "+String(intensity,12).fillLeft(' ',12)).c_str());
		}
		statusBar()->repaint();
	}
	
	void SpectrumMDIWindow::showGridLines(bool b)
	{
		SpectrumWindow* window = activeWindow();
		if (window!=0)
		{
			window->widget()->canvas()->showGridLines(b);
		}
	}
	
	void SpectrumMDIWindow::resetZoom()
	{
		SpectrumWindow* window = activeWindow();
		if (window!=0)
		{
			window->widget()->canvas()->resetZoom();
		}
	}
	
	void SpectrumMDIWindow::setActionMode(QAction* a)
	{
		SpectrumWindow* window = activeWindow();
		if (window!=0)
		{
			window->widget()->setActionMode(a);
		};
	}
	
	void SpectrumMDIWindow::setDrawMode(QAction* a)
	{
		Spectrum1DWindow* window = active1DWindow();
		if (window!=0)
		{
			window->widget()->canvas()->setDrawMode(a);
		}
	}
	
	void SpectrumMDIWindow::showPoints(bool on)
	{
		if (Spectrum2DWindow* win = active2DWindow())
		{
			win->widget()->canvas()->showPoints(on);
		}
	}
	
	void SpectrumMDIWindow::showColors(bool on)
	{
		if (Spectrum2DWindow* win = active2DWindow())
		{
			win->widget()->canvas()->showColors(on);
		}
	}
	
	void SpectrumMDIWindow::showContours(bool on)
	{
		if (Spectrum2DWindow* win = active2DWindow())
		{
			win->widget()->canvas()->showContours(on);
		}
	}
	
	void SpectrumMDIWindow::setIntensityScaledDots(bool on)
	{
		if (Spectrum2DWindow* win = active2DWindow())
		{
			win->widget()->canvas()->setIntensityScaledDots(on);
		}
	}

	void SpectrumMDIWindow::setBackView3D(bool on)
	{
		if(on)
		{
			show_back_view_3d_->setOn(false);	
			show_top_view_3d_->setOn(false);	
			intensity_scaled_dots_button_3d_->setOn(false);
			show_reset_view_3d_->setOn(false);
			if (Spectrum3DWindow* win = active3DWindow())
				{
					win->widget()->canvas()->openglwidget()->setBackView();
				}
		}
	}
	void SpectrumMDIWindow::setTopView3D(bool on)
	{
		if(on)
		{
			show_back_view_3d_->setOn(false);	
			show_top_view_3d_->setOn(false);	
			intensity_scaled_dots_button_3d_->setOn(false);
			show_reset_view_3d_->setOn(false);
			if (Spectrum3DWindow* win = active3DWindow())
			{
				win->widget()->canvas()->openglwidget()->setTopView();
			}
		}
	}
	void SpectrumMDIWindow::setIntensityScaledDots3D(bool on)
	{
			if(on)
			{
				show_back_view_3d_->setOn(false);	
				show_top_view_3d_->setOn(false);	
				intensity_scaled_dots_button_3d_->setOn(false);
				show_reset_view_3d_->setOn(false);
				
				if (Spectrum3DWindow* win = active3DWindow())
				{
						win->widget()->canvas()->openglwidget()->setIntensityScale(on);
				}
			}
	}

	void SpectrumMDIWindow::update3DToolbar(QWidget* w)
	{	
		if (Spectrum3DCanvas* wi = dynamic_cast<Spectrum3DCanvas*>(w))
		{
			set_pick_action_3d_->setOn(wi->openglwidget()->getShowSelect());
			set_zoom_action_3d_->setOn(wi->openglwidget()->getShowZoom());
		}
	}
	void SpectrumMDIWindow::update2DToolbar(QWidget* w)
	{
		if (Spectrum2DCanvas* wi = dynamic_cast<Spectrum2DCanvas*>(w))
		{
			show_points_button_2d_->setOn(wi->getShowPoints());
			show_colors_button_2d_->setOn(wi->getShowColors());
			show_contours_button_2d_->setOn(wi->getShowContours());
		}
	}
	
	void SpectrumMDIWindow::update1DToolbar(QWidget* w)
	{
		if (Spectrum1DCanvas* wi = dynamic_cast<Spectrum1DCanvas*>(w))
		{
			switch (wi->getDrawMode())
			{
				case Spectrum1DCanvas::DM_PEAKS:
					set_peak_mode_->setOn(true);
					break;
				case Spectrum1DCanvas::DM_CONNECTEDLINES:
					set_connected_lines_mode_->setOn(true);
					break;
				default:
					throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			};
		}
	}
	
	
	void SpectrumMDIWindow::updateToolbar(QWidget* w)
	{
		if (w)
		{
			//set draw mode
			if (dynamic_cast<Spectrum1DWindow*>(w)) 
			{
	
				switch (((Spectrum1DWindow*)w)->widget()->canvas()->getDrawMode())
				{
					case Spectrum1DCanvas::DM_PEAKS:
						set_peak_mode_->setOn(true);
						break;
					case Spectrum1DCanvas::DM_CONNECTEDLINES:
						set_connected_lines_mode_->setOn(true);
						break;
					default:
						throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				};
	
				//set action mode
				switch (((Spectrum1DWindow*)w)->widget()->getActionMode())
				{
					case SpectrumCanvas::AM_SELECT:
						set_pick_action_->setOn(true);
						break;
					case SpectrumCanvas::AM_ZOOM:
						set_zoom_action_->setOn(true);
						break;
					case SpectrumCanvas::AM_TRANSLATE:
						set_translate_action_->setOn(true);
						break;
					default:
						throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				};
	
				//set grid mode
				grid_button_->setOn(((Spectrum1DWindow*)w)->widget()->canvas()->gridLinesShown());
	
				//update link selector
				QFileInfo file;
				int current_index=0;
				int active_linked_to_address = link_map_[PointerSizeInt(&(*ws_->activeWindow()))];
				link_box_->clear();
				link_box_->insertItem("<unlinked>");
				QWidgetList windows = ws_->windowList();
				for ( int i = 0; i < int(windows.count()); ++i )
				{
					QWidget *window = windows.at(i);
					if (window!=w)
					{
						current_index++;
						file = QFileInfo(window->caption());
						link_box_->insertItem(file.fileName()+"  ("+window->caption()+")");
						if (active_linked_to_address==PointerSizeInt(&(*window)))
							{
								link_box_->setCurrentItem(current_index);
							}
					}
				}
				tool_bar_2d_->hide();
				tool_bar_3d_->hide();
				tool_bar_->show();
			}
			else if (dynamic_cast<Spectrum2DWindow*>(w))
			{
				Spectrum2DWindow* wi = dynamic_cast<Spectrum2DWindow*>(w);
				show_points_button_2d_->setOn(wi->widget()->canvas()->getShowPoints());
				show_colors_button_2d_->setOn(wi->widget()->canvas()->getShowColors());
				show_contours_button_2d_->setOn(wi->widget()->canvas()->getShowContours());
	      intensity_scaled_dots_button_2d_->setOn(wi->widget()->canvas()->isIntensityScaledDots());
	
				//set action mode
				switch (wi->widget()->getActionMode())
				{
					case SpectrumCanvas::AM_SELECT:
						set_pick_action_2d_->setOn(true);
						break;
					case SpectrumCanvas::AM_ZOOM:
						set_zoom_action_2d_->setOn(true);
						break;
					case SpectrumCanvas::AM_TRANSLATE:
						set_translate_action_2d_->setOn(true);
						break;
					case SpectrumCanvas::AM_MEASURE:
						set_measure_action_2d_->setOn(true);
						break;
					default:
						throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				};

				//set grid mode
				grid_button_2d_->setOn(((Spectrum2DWindow*)w)->widget()->canvas()->gridLinesShown());

				tool_bar_->hide();
				tool_bar_3d_->hide();
				tool_bar_2d_->show();
			}
			else if (dynamic_cast<Spectrum3DWindow*>(w))
			{
				Spectrum3DWindow* wi = dynamic_cast<Spectrum3DWindow*>(w);
				//set action mode
				switch (wi->widget()->getActionMode())
				{
					case SpectrumCanvas::AM_SELECT:
						set_pick_action_3d_->setOn(true);
						break;
					case SpectrumCanvas::AM_ZOOM:
						set_zoom_action_3d_->setOn(true);
						break;
					default:
						throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				};			
				tool_bar_2d_->hide();
				tool_bar_->hide();
				tool_bar_3d_->show();
			}
			//layer manager
			updateLayerbar();
		}
	}
	
	void SpectrumMDIWindow::updateLayerbar()
	{
	
		layer_bar_->hide();
		layer_manager_->disconnect(SIGNAL(visibilityChanged(int, bool)));
		layer_manager_->disconnect(SIGNAL(activatedChanged(int)));
		layer_manager_->disconnect(SIGNAL(removed(int)));
		layer_manager_->reset();

		SpectrumCanvas* cc;
		if (active1DWindow() != 0 )
		{
			cc = active1DWindow()->widget()->canvas();
	  }
	  else if (active2DWindow() != 0 )
		{
			cc = active2DWindow()->widget()->canvas();
		}
		else if (active3DWindow()!= 0 )
		{
			cc = active3DWindow()->widget()->canvas();
		}
		else
		{
			return;
		}
		
		for (UnsignedInt i = 0; i<cc->getDataSetCount(); ++i)
		{
			layer_manager_->addLayer(cc->getDataSetName(i));
			layer_manager_->setVisible(i,cc->isDataSetVisible(i));
	 	}
		
		layer_manager_->activate(cc->activeDataSetIndex());
	
		connect(layer_manager_,SIGNAL(visibilityChanged(int, bool)),cc,SLOT(changeVisibility(int, bool)));
		connect(layer_manager_,SIGNAL(activatedChanged(int)),cc,SLOT(activateDataSet(int)));
		connect(layer_manager_,SIGNAL(removed(int)),cc,SLOT(removeDataSet(int)));
		layer_bar_->show();
	}
	
	
	void SpectrumMDIWindow::updateTabBar(QWidget* w)
	{
		if (w)
		{
			tab_bar_->setCurrentTab(PointerSizeInt(&(*w)));
		}
	}
	
	void SpectrumMDIWindow::tileVertical()
	{
	    // primitive horizontal tiling
	    QWidgetList windows = ws_->windowList();
	    if ( !windows.count() )
	        return;
	
	    if (active1DWindow())
	  	  active1DWindow()->showNormal();
	    if (active2DWindow())
	  	  active2DWindow()->showNormal();
	
	    int heightForEach = ws_->height() / windows.count();
	    int y = 0;
	    for ( int i = 0; i < int(windows.count()); ++i )
	    {
	        QWidget *window = windows.at(i);
	        if ( window->isMaximized() || window->isFullScreen() )
	        {
	
	            // prevent flicker
	            window->hide();
	            window->setWindowState(Qt::WindowNoState);
	            window->show();
	        }
	        int preferredHeight = window->minimumHeight()+window->parentWidget()->baseSize().height();
	        int actHeight = QMAX(heightForEach, preferredHeight);
	
	        window->parentWidget()->setGeometry( 0, y, ws_->width(), actHeight );
	        y += actHeight;
	    }
	}
	
	void SpectrumMDIWindow::tileHorizontal()
	{
	    // primitive horizontal tiling
	    QWidgetList windows = ws_->windowList();
	    if ( !windows.count() )
	        return;
	
	    if (active1DWindow())
	  	  active1DWindow()->showNormal();
	    if (active2DWindow())
	  	  active2DWindow()->showNormal();
	
	    int widthForEach = ws_->width() / windows.count();
	    int y = 0;
	    for ( int i = 0; i < int(windows.count()); ++i )
	    {
	        QWidget *window = windows.at(i);
	        if ( window->testWState( WState_Maximized ) )
	        {
	            // prevent flicker
	            window->hide();
	            window->showNormal();
	        }
	        int preferredWidth = window->minimumWidth()+window->parentWidget()->baseSize().width();
	
	        int actWidth = QMAX(widthForEach, preferredWidth);
	
	        window->parentWidget()->setGeometry( y, 0, actWidth , ws_->height() );
	        y += actWidth;
	    }
	}
	
	void SpectrumMDIWindow::removeWidgetFromBar(QObject* o)
	{
		tab_bar_->removeTab(tab_bar_->tab(PointerSizeInt(&(*o))));
		id_map_.erase(PointerSizeInt(&(*o)));
	}
	
	void SpectrumMDIWindow::connectWindowSignals(SpectrumWindow* sw)
	
	{
		connect(sw,SIGNAL(sendStatusMessage(std::string,OpenMS::UnsignedInt)),this,SLOT(showStatusMessage(std::string,OpenMS::UnsignedInt)));
		connect(sw,SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
		connect(sw,SIGNAL(modesChanged(QWidget*)),this,SLOT(updateToolbar(QWidget*)));
		connect(sw,SIGNAL(openPreferences()),this,SLOT(preferencesDialog()));
		connect(sw,SIGNAL(destroyed()),this,SLOT(windowClosed()));
	}
	
	//! returns selected peaks of active spectrum framed by \c data_set_.begin() and the last peak BEFORE \c data_set_.end();
	vector<MSExperiment<>::SpectrumType::Iterator> SpectrumMDIWindow::getActiveSpectrumSelectedPeaks()
	{
		Spectrum1DWindow* w1 = active1DWindow();
		if (w1)
		{
			return (w1->widget()->canvas()->getSelectedPeaks());
		}
		return vector<MSExperiment<>::SpectrumType::Iterator>();
	}
	
	void SpectrumMDIWindow::gotoDialog()
	{
		SpectrumWindow* w = activeWindow();
		if (w)
		{
			w->showGoToDialog();
		}
	}
	
	
	void SpectrumMDIWindow::showPeaklistActiveSpectrum()
	{
		std::vector<MSExperiment<>::SpectrumType::Iterator> peaks = getActiveSpectrumSelectedPeaks();
	
		//only if peaks selected:
		if (peaks.size() >= 3)
		{
			std::ostringstream peaklist;
	
			for (std::vector<MSExperiment<>::SpectrumType::Iterator>::iterator it = peaks.begin(); it != peaks.end(); it++)
			{
				if (it == peaks.begin())
				{
					peaklist << "First Peak in Spectrum: Position " << (**it).getPosition()[0] << ", Intensity " << (**it).getIntensity() << endl;
				}
				else if (*it == *(peaks.rbegin()))
				{
					peaklist << "Last Peak in Spectrum: Position " << (**it).getPosition()[0] << ", Intensity " << (**it).getIntensity();
				}
				else
				{
					peaklist << "Selected Peak: Position " << (**it).getPosition()[0] << ", Intensity " << (**it).getIntensity() << endl;
				}
			}
			QMessageBox::information(this, QString("Peaklist"), QString(peaklist.str().c_str()));
		}
		else
		{
			QMessageBox::information(this, QString("Peaklist"), QString("No Peaks selected!"));
		}
	
	}
	
	void SpectrumMDIWindow::pickActiveSpectrum()
	{
		Spectrum1DWindow* w = active1DWindow();
		if (w!=0)
		{
			//handle intensity mods
			bool switched = false;
			if (w->widget()->canvas()->getIntensityMode() == SpectrumCanvas::IM_LOG)
			{
				w->widget()->setIntensityMode(SpectrumCanvas::IM_NONE);
				switched = true;
			}
			
			//pick data
			DPeakPickerCWT<1> peak_picker;
			MSExperiment<DPickedPeak<1> > picked;
			peak_picker(picked);
			MSExperiment<DRawDataPoint<1> > raw;
			raw.resize(1);
			raw[0].insert(raw[0].begin(), w->widget()->canvas()->currentDataSet()[0].begin(), w->widget()->canvas()->currentDataSet()[0].end());
			peak_picker.pick(raw);
	
			//add new spectrum
			String new_name = w->widget()->canvas()->currentDataSet().getName()+" (picked)";
			Spectrum1DCanvas::ExperimentType& exp = w->widget()->canvas()->addEmptyDataSet();
			exp.setName(new_name); // set layername
			exp.resize(1);
			exp[0].insert(exp[0].begin(), picked[0].begin(), picked[0].end());
			//color picked peaks
			for (Spectrum1DCanvas::ExperimentType::SpectrumType::Iterator it = exp[0].begin(); it!= exp[0].end(); ++it)
			{
				it->setMetaValue(UnsignedInt(5),string("#FF00FF"));
			}
			
			w->widget()->canvas()->finishAdding();
			
			updateLayerbar();
	
			if (switched)
			{
				w->widget()->setIntensityMode(SpectrumCanvas::IM_LOG);
			}
		}
	}
	
	SpectrumWindow*  SpectrumMDIWindow::activeWindow() const
	{
		return dynamic_cast<SpectrumWindow*>(ws_->activeWindow());
	}
	
	Spectrum1DWindow* SpectrumMDIWindow::active1DWindow() const
	{
		Spectrum1DWindow* s1;
		if (s1 = dynamic_cast<Spectrum1DWindow*>(ws_->activeWindow()))
		{
			return s1;
		}
		return 0;
	}
	
	Spectrum2DWindow* SpectrumMDIWindow::active2DWindow() const
	{
		Spectrum2DWindow* s2;
		if (s2 = dynamic_cast<Spectrum2DWindow*>(ws_->activeWindow()))
		{
			return s2;
		}
		return 0;
	}
	
	Spectrum3DWindow* SpectrumMDIWindow::active3DWindow() const
	{
		Spectrum3DWindow* s3;
		if (s3 = dynamic_cast<Spectrum3DWindow*>(ws_->activeWindow()))
		{
			return s3;
		}
		return 0;
	}
	
	PreferencesDialogPage* SpectrumMDIWindow::createPreferences(QWidget* parent)
	{
		return new SpectrumMDIWindowPDP(this,parent);
	}
	
	void SpectrumMDIWindow::loadPreferences(string filename)
	{
		//compose default ini file path
		String default_ini_file;
	  char * home;
	  home = getenv ("HOME");
	  if (home!=NULL)
	  {
	  	default_ini_file = home;
	  	default_ini_file = default_ini_file + "/";
	  }
		default_ini_file = default_ini_file + ".TOPPView.ini";
		
		if (filename=="")
		{
			filename = default_ini_file;
		}
		prefs_.setValue("PreferencesFile" , filename);
		
		//load preferences, if file exists
		FILE * infile;
	  infile = fopen (filename.c_str(), "r");
	  if (infile != NULL)
		{
			prefs_.clear();
			prefs_.setValue("PreferencesFile" , filename);
			prefs_.load(filename);
	  }
		else
		{
			if (filename != default_ini_file)
			{
				cerr << "Unable to load INI File: '" << filename << "'" << endl;
			}
			close (infile);
		}
		
		//set missing defaults
	  checkPreferences_();
	  
		//set the recent files
		Param p = prefs_.copy("Preferences:RecentFiles");
		if (p.size()!=0)
		{
			for (Param::ConstIterator it=p.begin() ; it!=p.end() ; ++it)
			{
				recent_files_.push_back(string(it->second));
			}
		}

		updateRecentMenu_();
	}
	
	void SpectrumMDIWindow::savePreferences()
	{
		// replace recent files
		prefs_.remove("Preferences:RecentFiles");
	
		for (UnsignedInt i=0; i<recent_files_.size(); ++i)
		{
			prefs_.setValue("Preferences:RecentFiles:"+String(i),recent_files_[i]);
		}
		
		//save only the subsection that begins with "Preferences:"
		try
		{
			prefs_.copy("Preferences:").store(string(prefs_.getValue("PreferencesFile")));
		}
		catch(Exception::UnableToCreateFile& e)
		{
			cerr << "Unable to create INI File: '" << string(prefs_.getValue("PreferencesFile")) << "'" << endl;
		}
	}
	
	
	void SpectrumMDIWindow::checkPreferences_()
	{
		Param default_preferences;
		
		//general
		default_preferences.setValue("DB:Host", "localhost");
		default_preferences.setValue("DB:Login", "NoName");
		default_preferences.setValue("DB:Name", "OpenMS");
		default_preferences.setValue("DB:Port", "3306");
		default_preferences.setValue("DefaultMapView", "2D");
		default_preferences.setValue("DefaultPath", ".");
		default_preferences.setValue("NumberOfRecentFiles", 15);
		default_preferences.setValue("Legend", "Show");	
		
		//1d
		default_preferences.setValue("1D:HighColor", "#ff0000");
		default_preferences.setValue("1D:IconColor", "#000000");
		default_preferences.setValue("1D:PeakColor", "#0000ff");
		default_preferences.setValue("1D:BackgroundColor", "#ffffff");
		default_preferences.setValue("1D:Mapping:MappingOfMzTo","X-Axis");
		
		//2d
		default_preferences.setValue("2D:BackgroundColor", "#ffffff");
		default_preferences.setValue("2D:MarchingSquaresSteps", 20);
		default_preferences.setValue("2D:InterpolationSteps", 200);
		default_preferences.setValue("2D:Dot:Gradient", "Linear|0,#ffff00;11,#ffaa00;32,#ff0000;55,#aa00ff;78,#5500ff;100,#000000");
		default_preferences.setValue("2D:Dot:Mode", 1);
		default_preferences.setValue("2D:Surface:Gradient", "Linear|0,#ffffff;22,#fdffcb;50,#ffb4b4;75,#d7cfff;100,#c1c1c1");
		default_preferences.setValue("2D:Contour:Lines", 8);
		default_preferences.setValue("2D:Mapping:MappingOfMzTo","X-Axis");
	
		//3d
		default_preferences.setValue("3D:Dot:Mode", 1);
		default_preferences.setValue("3D:Shade:Mode", 1);
		default_preferences.setValue("3D:Dot:Gradient", "Linear|0,#ffff00;11,#ffaa00;32,#ff0000;55,#aa00ff;78,#5500ff;100,#000000");
		default_preferences.setValue("3D:Dot:InterpolationSteps",200);
		default_preferences.setValue("3D:BackgroundColor", "#ffffff");
		default_preferences.setValue("3D:AxesColor", "#000000");
		default_preferences.setValue("3D:IntMode",0);	
		default_preferences.setValue("3D:Dot:LineWidth",2);
		default_preferences.setValue("3D:IntScale:Mode",0);
		prefs_.setDefaults(default_preferences,"Preferences");
	}
	

	void SpectrumMDIWindow::openRecentFile(int i)
	{
		if (getPrefAsString("Preferences:DefaultMapView")=="2D")
		{
			addSpectrum(recent_files_[i].c_str(),true,true,true);
		}
		else
		{
			addSpectrum(recent_files_[i].c_str(),true,false,true);
		}
	}
	
	void SpectrumMDIWindow::findFeaturesActiveSpectrum()
		{
			Spectrum2DWindow* w = active2DWindow();
			if (w!=0)
			{
				FeaFiDialog dialog(this, "FeaFiDialog");
				if (dialog.exec() == QDialog::Accepted)
				{
					//handle intensity mod
					bool switched = false;
					if (w->widget()->canvas()->getIntensityMode() == SpectrumCanvas::IM_LOG)
					{
						w->widget()->setIntensityMode(SpectrumCanvas::IM_NONE);
						switched = true;
					}
					
					//find features
					FeatureFinder& finder = dialog.getFeatureFinder();
					Spectrum2DCanvas::ExperimentType in = w->widget()->canvas()->currentDataSet();
					finder.setData(in);
					Spectrum2DCanvas::ExperimentType out;
					//copy to sort features by RT
					DFeatureMap<2> features = finder.run();
					features.sortByPosition();
					out.set2DData(features);
					
					//display features
					setFeatureMap_(w->widget()->canvas(), out, w->widget()->canvas()->currentDataSet().getName());
					updateLayerbar();
					
					//handle intensity mode
					if (switched)
					{
						w->widget()->setIntensityMode(SpectrumCanvas::IM_LOG);
					}
				}
			}
		}

	void SpectrumMDIWindow::setFeatureMap_(Spectrum2DCanvas* canvas, Spectrum2DCanvas::ExperimentType& exp, String caption)
	{
		if (exp.size()==0) 
		{
			return;
		}
	
		// First layer: just the centers
		exp.setName("Center (Features of " +caption+")");
		canvas->addDataSet(exp);
	
		// Second layer: convex hulls
		exp.setName("Convex hulls (Features of " +caption+ ")");
		exp.setMetaValue("FeatureDrawMode",std::string("ConvexHulls"));
		canvas->addDataSet(exp);
	}	
	
	void SpectrumMDIWindow::closeEvent(QCloseEvent * e)
	{
		savePreferences();
		e->accept();
	}
	
	void SpectrumMDIWindow::windowClosed()
	{
		// close tab bar, when last window is closed
		if (ws_->windowList().count()==1)
		{
			layer_bar_->hide();
		}
	}


	void SpectrumMDIWindow::openSpectrumDialog()
	{
		OpenDialog dialog(prefs_,this,"Open Spectrum Dialog");
		if (dialog.exec())
		{
			//Open Files
			if (dialog.getSource()==OpenDialog::FILE)
			{
				for(vector<String>::const_iterator it=dialog.getNames().begin();it!=dialog.getNames().end();it++)
				{
					addSpectrum(*it,dialog.isOpenAsNewTab(),dialog.isViewMaps2D(),true,dialog.getMower());
				}
			}
			else
			// Open from DB
			{
				for(vector<String>::const_iterator it=dialog.getNames().begin();it!=dialog.getNames().end();it++)
				{
					addDBSpectrum(it->toInt(),dialog.isOpenAsNewTab(),dialog.isViewMaps2D(),true,dialog.getMower());
				}
			}
			maximizeActiveSpectrum();
		}
	}

} //namespace OpenMS

