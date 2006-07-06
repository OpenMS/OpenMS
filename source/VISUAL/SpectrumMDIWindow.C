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
#include <OpenMS/FORMAT/FileHandler.h>
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

//action modes
#include "ICONS/zoom.xpm"
#include "ICONS/translate.xpm"
#include "ICONS/select.xpm"
#include "ICONS/measure.xpm"

//intensity modes
#include "ICONS/lin.xpm"
#include "ICONS/log.xpm"
#include "ICONS/percentage.xpm"
#include "ICONS/snap.xpm"

//common
#include "ICONS/grid.xpm"
#include "ICONS/print.xpm"
#include "ICONS/reset_zoom.xpm"
#include "ICONS/tile_horizontal.xpm"
#include "ICONS/tile_vertical.xpm"

//1d
#include "ICONS/lines.xpm"
#include "ICONS/peaks.xpm"

//2d
#include "ICONS/points.xpm"
#include "ICONS/colors.xpm"
#include "ICONS/contours.xpm"



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
		windows->insertItem(QIconSet(QPixmap(XPM_tile_h)),"Tile &vertical",this,SLOT(tileHorizontal()));
		windows->insertItem(QIconSet(QPixmap(XPM_tile_v)),"Tile &horizontal",this,SLOT(tileVertical()));
		//Tools menu
		tools_menu_ = new QPopupMenu(this);
		menuBar()->insertItem("&Tools", tools_menu_);
		tools_menu_->insertItem("&Show Selected Peaks (1D)", this, SLOT(showPeaklistActiveSpectrum()));
		tools_menu_->insertItem("&Pick Peaks (1D)", this, SLOT(pickActiveSpectrum()));
		tools_menu_->insertItem("&Find Features (2D)", this, SLOT(findFeaturesActiveSpectrum()));
	
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
	
		//create toolbars and connect signals
		createToolBar_();
	
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
	
	
	void SpectrumMDIWindow::addDBSpectrum(UnsignedInt db_id, bool as_new_window, bool maps_as_2d, bool maximize, OpenDialog::Mower use_mower)
	{
	#ifdef DB_DEF
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
		
		SpectrumCanvas::ExperimentType* exp;
		MSExperiment<>* exp2; //temporary data

		if (activeWindow_()==0)
		{
			as_new_window = true;
		}
		
		//open in new window
		if (as_new_window)
		{
			//create 1D View
			if (con.lastResult().value(0).toInt()==1)
			{
				// create 1D window
				w = new Spectrum1DWindow(ws_,"Spectrum1DWindow",WDestructiveClose);
				
				//determine Spectrum id
				con.executeQuery("SELECT id from Spectrum where MSExperiment_id='"+db_id_string+"' and MS_Level='1'");
				con.lastResult().first();
				UID spectrum_id = con.lastResult().value(0).toInt();
				
				//load data
				exp = &(w->widget()->canvas()->addEmptyDataSet());
				exp2->push_back(* (dba.loadSpectrum(spectrum_id)));
			}
			//create 2D/3D view
			else
			{
				//create 2D view
				if (maps_as_2d)
				{
					//create 2D window
					w = new Spectrum2DWindow(ws_,"Spectrum2DWindow",WDestructiveClose);
					
					//load spectrum
					exp2 = dba.loadMSExperiment(db_id);
					exp = &(w->widget()->canvas()->addEmptyDataSet());
				}
				//create 3D view
				else
				{
					// create 3D window
					w = new Spectrum3DWindow(ws_, "Spectrum3DWindow", WDestructiveClose);		
					
					//load data
					exp2 = dba.loadMSExperiment(db_id);
					exp = &(w->widget()->canvas()->addEmptyDataSet());
				}
				w->widget()->setMainPreferences(prefs_);
			}
		}
		//open in active window
		else
		{
			//create 1D View
			if (con.lastResult().value(0).toInt()==1)
			{
				w = active1DWindow_();
				//wrong active window type
				if (w==0)
				{
					QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+db_id_string+") in a 2D window!<BR>Please open the file in new tab.").c_str());
					return;
				}
				else //open it
				{
					//determine Spectrum id
					con.executeQuery("SELECT id from Spectrum where MSExperiment_id='"+db_id_string+"' and MS_Level='1'");
					con.lastResult().first();
					UID spectrum_id = con.lastResult().value(0).toInt();
					
					//load data
					exp = &(w->widget()->canvas()->addEmptyDataSet());
					exp2->push_back(*(dba.loadSpectrum(spectrum_id)));
				}
			}
			//create 2D/3D view
			else
			{
				Spectrum2DWindow* w2 = active2DWindow_();
				Spectrum3DWindow* w3 = active3DWindow_();
				//wrong active window type
				if (w2==0 && w3==0)
				{
					QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+db_id_string+") in a 2D/3D window!<BR>Please open the file in new tab.").c_str());
					return;
				}
				//create 2D view
				if (w2!=0)
				{
					w = w2;
					
					//load data
					MSExperiment<>* exp = dba.loadMSExperiment(db_id);
					exp = &(w->widget()->canvas()->addEmptyDataSet());
					*exp = *exp2;
					delete(exp2);
				}
				//create 3D view
				else
				{
					w = w3;
					
					//load data
					exp = dba.loadMSExperiment(db_id);
					exp = &(w->widget()->canvas()->addEmptyDataSet());
				}
			}

			*exp = *exp2;
			delete(exp2);
			
			//noise estimator
			float cutoff = 0;
			if(use_mower!=OpenDialog::NO_MOWER && exp->size()>1)
			{
			  cutoff = estimateNoise_(*exp);
			}
			
			exp->setName(caption);
			w->widget()->canvas()->finishAdding();

			//use_mower
			
			//do for all windows
			if (as_new_window)
			{
				w->widget()->canvas()->finishAdding(cutoff);
				connectWindowSignals_(w);
				w->setCaption(caption.c_str());
				addClient(w,caption);
				addTab_(w,caption);
			}
			
			//do for all (in active and in new window, 1D/2D/3D)
			if(maximize)
			{
				w->showMaximized();
			}

			//do for all windows
			updateLayerbar();
		}
	#endif
	}

	float SpectrumMDIWindow::estimateNoise_(const SpectrumCanvas::ExperimentType& exp)
	{
		float noise = 0.0;
		UnsignedInt count = 0;
		srand(time(0));
		while (count<10)
		{	
			UnsignedInt scan = (UnsignedInt)( (double)rand() / (double)(RAND_MAX+1) * (exp.size()-1) );
			 
			if (scan < exp.size() && exp[scan].getMSLevel()==1)
			{
				vector<float> tmp;
				tmp.reserve(exp[scan].size());
				for(SpectrumCanvas::ExperimentType::SpectrumType::ConstIterator it = exp[scan].begin()
					  ; it != exp[scan].end()
					  ; ++it)
				{
					tmp.push_back(it->getIntensity());
				}
				std::sort(tmp.begin(),tmp.end());
				noise += tmp[(UnsignedInt)ceil((float)(tmp.size()-1)/1.25)];
				
				++count;
			}
		}
		noise /= 10.0;
		return noise;
	}
	
	void SpectrumMDIWindow::preferencesDialog()
	{
		setActive(true);
		if (showPreferencesDialog())
		{
			savePreferences();
		}
	}

	void SpectrumMDIWindow::addSpectrum(const String& filename,bool as_new_window, bool maps_as_2d, bool maximize, OpenDialog::Mower use_mower, FileHandler::Type force_type)
	{
		//for comparison
		String filename_lower(filename);
		filename_lower.toLower();
	  	
		if (filename == "")
		{
			return;
		}
		QFileInfo file(filename.c_str());
		if (!file.exists())
		{
			QMessageBox::warning(this,"Open file error",("The file '"+filename+"' does not exist!").c_str());
			return;
		}

		//extract the filename without path
		String caption = file.fileName().ascii();
		
		//update recent files list
		addRecentFile_(filename);

		//windowpointer
		SpectrumWindow* w=0;
		
		if (activeWindow_()==0)
		{
			as_new_window = true;
		}
		
		if (as_new_window)
		{
			if (filename_lower.hasSuffix(".dta")) //1D
			{				
				w = new Spectrum1DWindow(ws_,"Spectrum1DWindow",WDestructiveClose);	
			}
			else if (maps_as_2d || filename_lower.hasSuffix(".feat")) //2d or features
			{
				w = new Spectrum2DWindow(ws_,"Spectrum1DWindow",WDestructiveClose);	
			}
			else //3d
			{
				w = new Spectrum3DWindow(ws_,"Spectrum1DWindow",WDestructiveClose);	
			}			

			//set main preferences
			w->widget()->setMainPreferences(prefs_);
		}
		else //!as_new_window
		{
			if (active1DWindow_()!=0) //1D window
			{
				if (!filename_lower.hasSuffix(".dta"))
				{
					QMessageBox::warning(this,"Wrong file type",("You cannot open 2D data ("+filename+") in a 1D window!").c_str());
					return;
				}

				w = active1DWindow_();
			}
			else if (active2DWindow_()!=0) //2d window
			{
				if (filename_lower.hasSuffix(".dta"))
				{
					QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+filename+") in a 2D window!").c_str());
					return;
				}
				w = active2DWindow_();
			}
			else if (active3DWindow_()!=0)//3d window
			{
				if (filename_lower.hasSuffix(".dta"))
				{
					QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+filename+") in a 3D window!").c_str());
					return;
				}
				w = active3DWindow_();
			}
		}
		
		//try to read the data from file (feature)
		if (filename_lower.hasSuffix(".feat"))
		{
			DFeatureMap<2> map;
			try
			{	
				DFeatureMapFile().load(filename,map);
			}
			catch(Exception::Base& e)
			{
				QMessageBox::warning(this,"Error",(String("Error while reading feature file: ")+e.what()).c_str());
				return;
			}
			map.setName(caption);
			w->widget()->canvas()->addDataSet(map);
		}
		else
		{
			//try to read the data from file (raw/peak data)
			SpectrumCanvas::ExperimentType* exp = &(w->widget()->canvas()->addEmptyDataSet());
			try
			{	
				FileHandler().loadExperiment(filename,*exp, force_type);
			}
			catch(Exception::Base& e)
			{
				QMessageBox::warning(this,"Error",(String("Error while reading data file: ")+e.what()).c_str());
				return;
			}
	
			//check if only one scan is in a 2d file
			if (as_new_window && active1DWindow_()==0 && exp->size()==1)
			{
				delete(w);
				w = new Spectrum1DWindow(ws_,"Spectrum1DWindow",WDestructiveClose);
				w->widget()->setMainPreferences(prefs_);
				exp = &(w->widget()->canvas()->addEmptyDataSet());
				FileHandler().loadExperiment(filename,*exp, force_type);
			}
	
			//do for all (in active and in new window, 1D/2D/3D)
			float cutoff=0;
			
			//calculate noise
			if(use_mower!=OpenDialog::NO_MOWER && exp->size()>1)
			{
			  cutoff = estimateNoise_(*exp);
			}
			exp->setName(caption);  // set layername
			w->widget()->canvas()->finishAdding(cutoff);
		}
		
		updateLayerbar();

		if(maximize)
		{
			w->showMaximized();
		}
		
		if (as_new_window)
		{
			connectWindowSignals_(w);
			w->setCaption(filename.c_str());
			addClient(w,filename);
			addTab_(w,caption);
		}
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
	
	
	void SpectrumMDIWindow::addTab_(SpectrumWindow* w, const String& tabCaption)
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
		SpectrumWindow* window = activeWindow_();
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
		SpectrumWindow* window = activeWindow_();
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
		
		//action modes
		action_modes_ = new QActionGroup(tool_bar_);
		action_modes_->setExclusive(TRUE);
		am_zoom_ = new QAction( QString("Action: Zoom"), QPixmap(XPM_zoom), NULL, Key_Z, action_modes_,"AM_ZOOM",TRUE);
		am_zoom_->addTo(tool_bar_);
		am_translate_ = new QAction( QString("Action: Translate"), QPixmap(XPM_translate), NULL, Key_T, action_modes_,"AM_TRANSLATE",TRUE);
		am_translate_->addTo(tool_bar_);
		am_select_ = new QAction( QString("Action: Select"), QPixmap(XPM_select), NULL, Key_S, action_modes_,"AM_SELECT",TRUE);
		am_select_->addTo(tool_bar_);
		am_measure_ = new QAction( QString("Action: Measure"), QPixmap(XPM_measure), NULL, Key_M, action_modes_,"AM_MEASURE",TRUE);
		am_measure_->addTo(tool_bar_);
		connect(action_modes_,SIGNAL(selected(QAction*)),this,SLOT(setActionMode(QAction*)));
		
		tool_bar_->addSeparator();
		
		//intensity modes
		intensity_modes_ = new QActionGroup(tool_bar_);
		intensity_modes_->setExclusive(TRUE);
		im_none_ = new QAction( QString("Intensity: Normal"), QPixmap(XPM_lin), NULL, Key_N, intensity_modes_,"IM_NONE",TRUE);
		im_none_->addTo(tool_bar_);
		im_log_ = new QAction( QString("Intensity: Logarithmic"), QPixmap(XPM_log), NULL, Key_L, intensity_modes_,"IM_LOG",TRUE);
		im_log_->addTo(tool_bar_);
		im_percentage_ = new QAction( QString("Intensity: Percentage"), QPixmap(XPM_percentage), NULL, Key_P, intensity_modes_,"IM_PERCENTAGE",TRUE);
		im_percentage_->addTo(tool_bar_);
		im_snap_ = new QAction( QString("Intensity: Snap to maximum displayed intensity"), QPixmap(XPM_snap), NULL, Key_A, intensity_modes_,"IM_SNAP",TRUE);
		im_snap_->addTo(tool_bar_);
		connect(intensity_modes_,SIGNAL(selected(QAction*)),this,SLOT(setIntensityMode(QAction*)));
		
		tool_bar_->addSeparator();
					
		//common buttons
		reset_zoom_button_ = new QToolButton(QIconSet(QPixmap(XPM_reset_zoom)),"Reset Zoom (Backspace)", "Reset Zoom", NULL, NULL , tool_bar_, "resetZoomButton");
		reset_zoom_button_->setAccel(Key_Backspace);
		connect(reset_zoom_button_,SIGNAL(clicked()),this,SLOT(resetZoom()));

		tool_bar_->addSeparator();

		grid_button_ = new QToolButton(QIconSet(QPixmap(XPM_grid)),"Show grid (G)","Show grid",NULL, NULL ,tool_bar_,"gridButton");
		grid_button_->setAccel(Key_G);
		grid_button_->setToggleButton(true);
		connect(grid_button_,SIGNAL(toggled(bool)),this,SLOT(showGridLines(bool)));		
		
		tool_bar_->addSeparator();
	
		print_button_ = new QToolButton(QIconSet(QPixmap(XPM_print)),"Print (Ctrl + P)","print",NULL, NULL,tool_bar_,"printButton");
		print_button_->setAccel(CTRL + Key_P);
		connect(print_button_,SIGNAL(clicked()),this,SLOT(print()));

		tool_bar_->resize(tool_bar_->sizeHint());
		tool_bar_->show();

		// 1d toolbar
		tool_bar_1d_ = new QToolBar(this,"toolbar");
		
		draw_modes_ = new QActionGroup(tool_bar_1d_);
		draw_modes_->setExclusive(TRUE);
		dm_peaks_1d_ = new QAction( QString("Stick mode"), QPixmap(XPM_peaks), NULL, Key_I, draw_modes_,"DM_PEAKS",TRUE);
		dm_peaks_1d_->addTo(tool_bar_1d_);
		dm_rawdata_1d_ = new QAction( QString("Raw data mode"), QPixmap(XPM_lines), NULL, Key_R, draw_modes_,"DM_CONNECTEDLINES",TRUE);
		dm_rawdata_1d_->addTo(tool_bar_1d_);
		connect(draw_modes_,SIGNAL(selected(QAction*)),this,SLOT(setDrawMode1D(QAction*)));	
		
		tool_bar_1d_->addSeparator();
		
		link_box_ = new QComboBox(tool_bar_1d_);
		QToolTip::add(link_box_,"Use this combobox to link two spectra.\nLinked spectra zoom in/out together");
		connect(link_box_,SIGNAL(activated(const QString&)),this,SLOT(linkActiveTo(const QString&)));
		
		// 2d toolbar
		tool_bar_2d_ = new QToolBar(this,"toolbar");
		
		dm_points_2d_ = new QToolButton(QIconSet(QPixmap(XPM_points)), "Show dots (D)", "Show dots", 0, 0, tool_bar_2d_, "showPoints");
		dm_points_2d_->setAccel(Key_D);
		dm_points_2d_->setToggleButton(true);
		connect(dm_points_2d_, SIGNAL(toggled(bool)), this, SLOT(showPoints(bool)));
	
		dm_surface_2d_ = new QToolButton(QIconSet(QPixmap(XPM_colors)), "Show colored surface (U)", "Show colored surface", 0, 0, tool_bar_2d_, "showSurface");
		dm_surface_2d_->setAccel(Key_U);
		dm_surface_2d_->setToggleButton(true);
		connect(dm_surface_2d_, SIGNAL(toggled(bool)), this, SLOT(showSurface(bool)));
	
		dm_contours_2d_ = new QToolButton(QIconSet(QPixmap(XPM_contours)), "Show contour lines (C)", "Show contour lines", 0, 0, tool_bar_2d_, "showContours");
		dm_contours_2d_->setAccel(Key_C);
		dm_contours_2d_->setToggleButton(true);
		connect(dm_contours_2d_, SIGNAL(toggled(bool)), this, SLOT(showContours(bool)));
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
					connect(activeWindow_()->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),dynamic_cast<SpectrumWindow*>(window)->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
					connect(dynamic_cast<SpectrumWindow*>(window)->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),activeWindow_()->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
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
			disconnect(id_map_[active_linked_to_address]->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)), activeWindow_()->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
			disconnect(activeWindow_()->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)), id_map_[active_linked_to_address]->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
			//remove from the map
			link_map_.erase(active_address);
			link_map_.erase(active_linked_to_address);
		}
	}
	
	void SpectrumMDIWindow::showStatusMessage(string msg,UnsignedInt time)
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
		SpectrumWindow* window = activeWindow_();
		if (window!=0)
		{
			window->widget()->canvas()->showGridLines(b);
		}
	}
	
	void SpectrumMDIWindow::resetZoom()
	{
		SpectrumWindow* window = activeWindow_();
		if (window!=0)
		{
			window->widget()->canvas()->resetZoom();
		}
	}
	
	void SpectrumMDIWindow::setActionMode(QAction* a)
	{
		SpectrumWindow* w = activeWindow_();
		if (w)
		{
			string name = a->name();
			if(name=="AM_SELECT") w->widget()->setActionMode(SpectrumCanvas::SpectrumCanvas::AM_SELECT); 
			else if (name=="AM_ZOOM") w->widget()->setActionMode(SpectrumCanvas::SpectrumCanvas::AM_ZOOM); 
			else if(name=="AM_TRANSLATE") w->widget()->setActionMode(SpectrumCanvas::SpectrumCanvas::AM_TRANSLATE); 
			else if(name=="AM_MEASURE") w->widget()->setActionMode(SpectrumCanvas::SpectrumCanvas::AM_MEASURE); 
			else throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		};
	}

	void SpectrumMDIWindow::setIntensityMode(QAction* a)
	{
		SpectrumWindow* w = activeWindow_();
		if (w)
		{
			string name = a->name();
			if(name=="IM_NONE") w->widget()->setIntensityMode(SpectrumCanvas::SpectrumCanvas::IM_NONE); 
			else if (name=="IM_LOG") w->widget()->setIntensityMode(SpectrumCanvas::SpectrumCanvas::IM_LOG); 
			else if(name=="IM_PERCENTAGE") w->widget()->setIntensityMode(SpectrumCanvas::SpectrumCanvas::IM_PERCENTAGE); 
			else if(name=="IM_SNAP") w->widget()->setIntensityMode(SpectrumCanvas::SpectrumCanvas::IM_SNAP); 
			else throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		};
	}
	
	void SpectrumMDIWindow::setDrawMode1D(QAction* a)
	{
		Spectrum1DWindow* w = active1DWindow_();
		if (w)
		{
			string name = a->name();
			if (name == "DM_PEAKS") w->widget()->canvas()->setDrawMode(Spectrum1DCanvas::DM_PEAKS);
			else if (name == "DM_CONNECTEDLINES") w->widget()->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
			else throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
	}
	
	void SpectrumMDIWindow::showPoints(bool on)
	{
		if (Spectrum2DWindow* win = active2DWindow_())
		{
			win->widget()->canvas()->showPoints(on);
		}
	}
	
	void SpectrumMDIWindow::showSurface(bool on)
	{
		if (Spectrum2DWindow* win = active2DWindow_())
		{
			win->widget()->canvas()->showSurface(on);
		}
	}
	
	void SpectrumMDIWindow::showContours(bool on)
	{
		if (Spectrum2DWindow* win = active2DWindow_())
		{
			win->widget()->canvas()->showContours(on);
		}
	}
	
	void SpectrumMDIWindow::updateToolbar(QWidget* /*widget*/)
	{
		SpectrumWindow* w = activeWindow_();
		
		if (w)
		{
			//set action mode
			switch (w->widget()->getActionMode())
			{
				case SpectrumCanvas::AM_SELECT:
					am_select_->setOn(true);
					break;
				case SpectrumCanvas::AM_ZOOM:
					am_zoom_->setOn(true);
					break;
				case SpectrumCanvas::AM_TRANSLATE:
					am_translate_->setOn(true);
					break;
				case SpectrumCanvas::AM_MEASURE:
					am_measure_->setOn(true);
					break;
				default:
					throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
			
			//set intensity mode
			switch (w->widget()->canvas()->getIntensityMode())
			{
				case SpectrumCanvas::IM_NONE:
					im_none_->setOn(true);
					break;
				case SpectrumCanvas::IM_LOG:
					im_log_->setOn(true);
					break;
				case SpectrumCanvas::IM_PERCENTAGE:
					im_percentage_->setOn(true);
					break;
				case SpectrumCanvas::IM_SNAP:
					im_snap_->setOn(true);
					break;
				default:
					throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
	
			//grid lines
			grid_button_->setOn(w->widget()->canvas()->gridLinesShown());
		}
		
		//1D
		Spectrum1DWindow* w1 = active1DWindow_();
		if (w1) 
		{
			//draw mode
			switch (w1->widget()->canvas()->getDrawMode())
			{
				case Spectrum1DCanvas::DM_PEAKS:
					dm_peaks_1d_->setOn(true);
					break;
				case Spectrum1DCanvas::DM_CONNECTEDLINES:
					dm_rawdata_1d_->setOn(true);
					break;
				default:
					throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			};

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
			
			//show/hide toolbars and buttons
			tool_bar_1d_->show();
			tool_bar_2d_->hide();
			am_measure_->setEnabled(false);
			am_select_->setEnabled(true);
		}
	
		//2d
		Spectrum2DWindow* w2 = active2DWindow_();
		if (w2) 
		{
			//draw modes
			dm_points_2d_->setOn(w2->widget()->canvas()->dotsAreShown());
			dm_surface_2d_->setOn(w2->widget()->canvas()->surfaceIsShown());
			dm_contours_2d_->setOn(w2->widget()->canvas()->contoursAreShown());

			//show/hide toolbars and buttons
			tool_bar_1d_->hide();
			tool_bar_2d_->show();
			am_measure_->setEnabled(true);
			am_select_->setEnabled(true);
		}

		//1D
		Spectrum3DWindow* w3 = active3DWindow_();
		if (w3) 
		{
			//show/hide toolbars and buttons
			tool_bar_1d_->hide();
			tool_bar_2d_->hide();
			am_measure_->setEnabled(false);
			am_select_->setEnabled(false);
		}
		
		//layer manager
		updateLayerbar();
	}
	
	void SpectrumMDIWindow::updateLayerbar()
	{
	
		layer_bar_->hide();
		layer_manager_->disconnect(SIGNAL(visibilityChanged(int, bool)));
		layer_manager_->disconnect(SIGNAL(activatedChanged(int)));
		layer_manager_->disconnect(SIGNAL(removed(int)));
		layer_manager_->reset();

		SpectrumCanvas* cc;
		if (active1DWindow_() != 0 )
		{
			cc = active1DWindow_()->widget()->canvas();
	  }
	  else if (active2DWindow_() != 0 )
		{
			cc = active2DWindow_()->widget()->canvas();
		}
		else if (active3DWindow_()!= 0 )
		{
			cc = active3DWindow_()->widget()->canvas();
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
	
	    if (active1DWindow_())
	  	  active1DWindow_()->showNormal();
	    if (active2DWindow_())
	  	  active2DWindow_()->showNormal();
	
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
	
	    if (active1DWindow_())
	  	  active1DWindow_()->showNormal();
	    if (active2DWindow_())
	  	  active2DWindow_()->showNormal();
	
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
	
	void SpectrumMDIWindow::connectWindowSignals_(SpectrumWindow* sw)
	{
		connect(sw->widget()->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(updateToolbar(QWidget*)));
		connect(sw,SIGNAL(sendStatusMessage(std::string,OpenMS::UnsignedInt)),this,SLOT(showStatusMessage(std::string,OpenMS::UnsignedInt)));
		connect(sw,SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
		connect(sw,SIGNAL(modesChanged(QWidget*)),this,SLOT(updateToolbar(QWidget*)));
		connect(sw,SIGNAL(openPreferences()),this,SLOT(preferencesDialog()));
		connect(sw,SIGNAL(destroyed()),this,SLOT(windowClosed()));
	}
	
	//! returns selected peaks of active spectrum framed by \c data_set_.begin() and the last peak BEFORE \c data_set_.end();
	vector<MSExperiment<>::SpectrumType::Iterator> SpectrumMDIWindow::getActiveSpectrumSelectedPeaks()
	{
		Spectrum1DWindow* w1 = active1DWindow_();
		if (w1)
		{
			return (w1->widget()->canvas()->getSelectedPeaks());
		}
		return vector<MSExperiment<>::SpectrumType::Iterator>();
	}
	
	void SpectrumMDIWindow::gotoDialog()
	{
		SpectrumWindow* w = activeWindow_();
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
		Spectrum1DWindow* w = active1DWindow_();
		if (w!=0)
		{
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
		}
	}
	
	SpectrumWindow*  SpectrumMDIWindow::activeWindow_() const
	{
		return dynamic_cast<SpectrumWindow*>(ws_->activeWindow());
	}
	
	Spectrum1DWindow* SpectrumMDIWindow::active1DWindow_() const
	{
		Spectrum1DWindow* s1;
		if ((s1 = dynamic_cast<Spectrum1DWindow*>(ws_->activeWindow())))
		{
			return s1;
		}
		return 0;
	}
	
	Spectrum2DWindow* SpectrumMDIWindow::active2DWindow_() const
	{
		Spectrum2DWindow* s2;
		if ((s2 = dynamic_cast<Spectrum2DWindow*>(ws_->activeWindow())))
		{
			return s2;
		}
		return 0;
	}
	
	Spectrum3DWindow* SpectrumMDIWindow::active3DWindow_() const
	{
		Spectrum3DWindow* s3;
		if ((s3 = dynamic_cast<Spectrum3DWindow*>(ws_->activeWindow())))
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
		default_preferences.setValue("MapIntensityCutoff", "None");	
		
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
		setCursor(Qt::WaitCursor);
		OpenDialog::Mower mow = OpenDialog::NO_MOWER;
		if ( getPrefAsString("Preferences:MapIntensityCutoff")=="Noise Estimator")
		{
			mow = OpenDialog::NOISE_ESTIMATOR;
		}
		addSpectrum(recent_files_[i].c_str(),true,getPrefAsString("Preferences:DefaultMapView")=="2D",true,mow);
		setCursor(Qt::ArrowCursor);
	}
	
	void SpectrumMDIWindow::findFeaturesActiveSpectrum()
	{
		Spectrum2DWindow* w = active2DWindow_();
		if (w!=0)
		{
			FeaFiDialog dialog(this, "FeaFiDialog");
			if (dialog.exec() == QDialog::Accepted)
			{
				//find features
				FeatureFinder& finder = dialog.getFeatureFinder();
				
				SpectrumCanvas::ExperimentType in = w->widget()->canvas()->currentDataSet();
				finder.setData(in);
				
				DFeatureMap<2> map = finder.run();
				
				//display features
				map.setName("Features: "+in.getName());
				w->widget()->canvas()->addDataSet(map);
				updateLayerbar();
			}
		}
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
			setCursor(Qt::WaitCursor);
			if (dialog.getSource()==OpenDialog::FILE)
			{
				for(vector<String>::const_iterator it=dialog.getNames().begin();it!=dialog.getNames().end();it++)
				{
					addSpectrum(*it,dialog.isOpenAsNewTab(),dialog.isViewMaps2D(),true,dialog.getMower(),dialog.forcedFileType());
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
			setCursor(Qt::ArrowCursor);
			maximizeActiveSpectrum();
		}
	}

} //namespace OpenMS

