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

#include <OpenMS/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DB/DBAdapter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/VISUAL/DIALOGS/SaveImageDialog.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum3DWidget.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolaySVDFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/BASELINE/TopHatFilter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>
#include <OpenMS/FORMAT/FeaturePairsXMLFile.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewPrefDialog.h>
#include <OpenMS/ANALYSIS/ID/IDFeatureMapper.h>
#include <OpenMS/ANALYSIS/ID/IDSpectrumMapper.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

//Qt
#include <QtGui/QToolBar>
#include <QtGui/QDockWidget>
#include <QtGui/QListWidget>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPrinter>
#include <QtGui/QStatusBar>
#include <QtGui/QToolButton>
#include <QtGui/QMessageBox>
#include <QtGui/QListWidgetItem>
#include <QtGui/QToolTip>
#include <QtGui/QFileDialog>
#include <QtGui/QPainter>
#include <QtGui/QPrintDialog>
#include <QtCore/QDir>
#include <QtCore/QDate>
#include <QtCore/QProcess>
#include <QtGui/QWhatsThis>
#include <QtGui/QInputDialog>
#include <QtGui/QTextEdit>

//action modes
#include "../VISUAL/ICONS/zoom.xpm"
#include "../VISUAL/ICONS/selection.xpm"

//intensity modes
#include "../VISUAL/ICONS/lin.xpm"
#include "../VISUAL/ICONS/logarithm.xpm"
#include "../VISUAL/ICONS/percentage.xpm"
#include "../VISUAL/ICONS/snap.xpm"

//common
#include "../VISUAL/ICONS/reset_zoom.xpm"
#include "../VISUAL/ICONS/tile_horizontal.xpm"
#include "../VISUAL/ICONS/tile_vertical.xpm"

//1d
#include "../VISUAL/ICONS/lines.xpm"
#include "../VISUAL/ICONS/peaks.xpm"

//2d
#include "../VISUAL/ICONS/precursors.xpm"
#include "../VISUAL/ICONS/projections.xpm"
#include "../VISUAL/ICONS/convexhulls.xpm"
#include "../VISUAL/ICONS/numbers.xpm"

//misc
#include "../VISUAL/ICONS/TOPPView.xpm"
#include "../VISUAL/ICONS/Oesterberg.xpm"

#include <algorithm>

using namespace std;

namespace OpenMS
{
  using namespace Internal;
	using namespace Math;

  TOPPViewBase::TOPPViewBase(QWidget* parent):
      QMainWindow(parent),
      DefaultParamHandler("TOPPViewBase")
  {
  	setWindowTitle("TOPPView");
    setWindowIcon(QIcon(toppview));
    //prevents errors caused by too small width,height values
    setMinimumSize(400,400);

    // create dummy widget (to be able to have a layout), Tab bar and workspace
    QWidget* dummy = new QWidget(this);
    setCentralWidget(dummy);
    QVBoxLayout* box_layout = new QVBoxLayout(dummy);
    tab_bar_ = new EnhancedTabBar(dummy);
    tab_bar_->setWhatsThis("Double-click tab to close it.");
    tab_bar_->addTab("dummy");
    tab_bar_->setMinimumSize(tab_bar_->sizeHint());
    tab_bar_->removeTab(0);

    //connect slots and sigals for selecting spectra
    connect(tab_bar_,SIGNAL(currentChanged(int)),this,SLOT(focusByTab(int)));
    connect(tab_bar_,SIGNAL(doubleClicked(int)),this,SLOT(closeByTab(int)));

    box_layout->addWidget(tab_bar_);
    ws_=new QWorkspace(dummy);
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateToolbar()));
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateTabBar(QWidget*)));
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateLayerbar()));
    box_layout->addWidget(ws_);

	//################## MENUS #################
    // File menu
    QMenu* file = new QMenu("&File",this);
    menuBar()->addMenu(file);
    file->addAction("&Open",this,SLOT(openSpectrumDialog()), Qt::CTRL+Qt::Key_O);
    file->addAction("&Close",this,SLOT(closeFile()));
    file->addSeparator();
    file->addAction("&Edit INI file", this, SLOT(editParamDialog()));
		file->addSeparator();
    
    QMenu* recent_menu = new QMenu("&Recent files", this);
    recent_as_new_layer_ = recent_menu->addAction("open as new layer");
    recent_as_new_layer_->setCheckable(true);
    recent_menu->addSeparator();
    //create the max mumber of recent files actions
  	recent_actions_.resize(20);
		for (UInt i = 0; i<20; ++i)
		{
			recent_actions_[i] = recent_menu->addAction("",this,SLOT(openRecentFile()));
			recent_actions_[i]->setVisible(false);
		}
  	file->addMenu(recent_menu);

    file->addSeparator();
    file->addAction("&Preferences",this, SLOT(preferencesDialog()));
    file->addAction("&Quit",qApp,SLOT(quit()), Qt::CTRL+Qt::Key_Q);
    
    //Layer menu
    QMenu* layer = new QMenu("&Layer",this);
    menuBar()->addMenu(layer);
    layer->addAction("&Save visible data",this,SLOT(saveLayer()), Qt::CTRL+Qt::Key_S);
    layer->addAction("&Rename layer",this,SLOT(renameLayer()), Qt::CTRL+Qt::Key_R);
    layer->addAction("&Edit metadata",this,SLOT(editMetadata()));
    layer->addAction("&Intensity distribution",this,SLOT(layerIntensityDistribution()));
		layer->addSeparator();
    layer->addAction("Apply &TOPP tool", this, SLOT(showTOPPDialog()), Qt::CTRL+Qt::Key_T);
    layer->addAction("&Annotate with identifiction", this, SLOT(annotateWithID()), Qt::CTRL+Qt::Key_A);
		layer->addSeparator();
    layer->addAction("&Preferences",this, SLOT(layerPreferencesDialog()));
    
    //View menu
    QMenu * view = new QMenu("&View",this);
    menuBar()->addMenu(view);
    view->addAction("&Go to",this,SLOT(gotoDialog()), Qt::CTRL+Qt::Key_G);
   	view->addAction("Show/Hide &axis legends",this,SLOT(changeAxisVisibility()));
   	view->addAction("Show/Hide &grid lines",this,SLOT(changeGridLines()));
   	
    //Image menu
    QMenu * image = new QMenu("&Image",this);
    menuBar()->addMenu(image);
    image->addAction("&Save to file",this,SLOT(saveImage()));
    image->addAction("&Print",this,SLOT(print()), Qt::CTRL+Qt::Key_P);

    //Windows menu
    QMenu * windows = new QMenu("&Windows", this);
    menuBar()->addMenu(windows);
    windows->addAction("&Cascade",this->ws_,SLOT(cascade()));
    windows->addAction("&Tile automatic",this->ws_,SLOT(tile()));
    windows->addAction(QIcon(QPixmap(tile_h)),"Tile &vertical",this,SLOT(tileHorizontal()));
    windows->addAction(QIcon(QPixmap(tile_v)),"Tile &horizontal",this,SLOT(tileVertical()));
	windows->addSeparator();
	
		//Help menu
		QMenu* help = new QMenu("&Help", this);
		menuBar()->addMenu(help);
		help->addAction(QWhatsThis::createAction(help));
		help->addSeparator();
		help->addAction("&About",this,SLOT(showAboutDialog()));
		
    //create status bar
    message_label_ = new QLabel(statusBar());
    statusBar()->addWidget(message_label_,1);

    rt_label_ = new QLabel("RT: 12345678", statusBar());
    rt_label_->setMinimumSize(rt_label_->sizeHint());
    rt_label_->setText("");
    statusBar()->addPermanentWidget(rt_label_,0);
    mz_label_ = new QLabel("m/z: 12345678", statusBar());
    mz_label_->setMinimumSize(mz_label_->sizeHint());
    mz_label_->setText("");
    statusBar()->addPermanentWidget(mz_label_,0);
    int_label_ = new QLabel("Int: 123456789012", statusBar());
    int_label_->setMinimumSize(int_label_->sizeHint());
    int_label_->setText("");
    statusBar()->addPermanentWidget(int_label_,0);

		//################## TOOLBARS #################
    //create toolbars and connect signals
  	QToolButton* b;
  	
  	//--Basic tool bar for all views--
    tool_bar_ = addToolBar("Basic tool bar");
    
    //action modes
    action_group_ = new QButtonGroup(tool_bar_);
    action_group_->setExclusive(true);
    
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(zoom));
    b->setToolTip("Action: Zoom + Translate");
    b->setShortcut(Qt::Key_Z);
    b->setCheckable(true);
    b->setWhatsThis("Action mode: Zoom + Translate<BR><BR>This mode allows to navigate in the data."
    								" The default is to zoom, Press the CTRL key for translation mode.<BR><BR>"
    								"A double-click resets the zoom.");
    action_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::AM_ZOOM);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(selection));
    b->setToolTip("Action: Select + Measure");
    b->setShortcut(Qt::Key_S);
    b->setCheckable(true);
    b->setWhatsThis("Action mode: Select + Measure<BR><BR>This mode allows to select peaks and"
    								" measure distances between peaks. The default is to select peaks. Press the"
    								" CTRL key for measurment mode.");
    action_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::AM_SELECT);
		tool_bar_->addWidget(b);

    connect(action_group_,SIGNAL(buttonClicked(int)),this,SLOT(setActionMode(int)));
    tool_bar_->addSeparator();
	   
    //intensity modes
    intensity_group_ = new QButtonGroup(tool_bar_);
    intensity_group_->setExclusive(true);
    
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(lin));
    b->setToolTip("Intensity: Normal");
    b->setShortcut(Qt::Key_N);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Normal<BR><BR>Intensity is displayed unmodified.");
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_NONE);
		tool_bar_->addWidget(b);
    
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(logarithm));
    b->setToolTip("Intensity: Logarithmic");
    b->setShortcut(Qt::Key_L);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Logarithmic<BR><BR>Intensity is displayed in a logarithmic scale.");
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_LOG);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(percentage));
    b->setToolTip("Intensity: Percentage");
    b->setShortcut(Qt::Key_P);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Percentage<BR><BR>Intensity is displayed as a percentage of the layer"
    								" maximum intensity. If only one layer is displayed this mode behaves like the"
    								" normal mode. If more than one layer is displayed intensities are aligned.");
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_PERCENTAGE);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(snap));
    b->setToolTip("Intensity: Snap to maximum displayed intensity");
    b->setShortcut(Qt::Key_A);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Snap to maximum displayed intensity<BR><BR> In this mode the"
    								" color gradient is adapted to the maximum currently displayed intensity.");
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_SNAP);
		tool_bar_->addWidget(b);

    connect(intensity_group_,SIGNAL(buttonClicked(int)),this,SLOT(setIntensityMode(int)));
    tool_bar_->addSeparator();

    //common buttons
    QAction* reset_zoom_button = tool_bar_->addAction(QPixmap(reset_zoom), "Reset Zoom", this, SLOT(resetZoom()));
    reset_zoom_button->setWhatsThis("Reset zoom: Zooms out as far as possible.");
    reset_zoom_button->setShortcut(Qt::Key_Backspace);

    tool_bar_->show();
    
    //--1D toolbar--
    tool_bar_1d_ = addToolBar("1D tool bar");

    //draw modes 1D
    draw_group_1d_ = new QButtonGroup(tool_bar_1d_);
    draw_group_1d_->setExclusive(true);
    
    b = new QToolButton(tool_bar_1d_);
    b->setIcon(QPixmap(peaks));
    b->setToolTip("Peak mode");
    b->setShortcut(Qt::Key_I);
    b->setCheckable(true);
    b->setWhatsThis("1D Draw mode: Peaks<BR><BR>Peaks are diplayed as sticks.");
    draw_group_1d_->addButton(b,Spectrum1DCanvas::DM_PEAKS);
		tool_bar_1d_->addWidget(b);
    
    b = new QToolButton(tool_bar_1d_);
    b->setIcon(QPixmap(lines));
    b->setToolTip("Raw data mode");
    b->setShortcut(Qt::Key_R);
    b->setCheckable(true);
    b->setWhatsThis("1D Draw mode: Raw data<BR><BR>Peaks are diplayed as a continous line.");
    draw_group_1d_->addButton(b,Spectrum1DCanvas::DM_CONNECTEDLINES);
		tool_bar_1d_->addWidget(b);

    connect(draw_group_1d_,SIGNAL(buttonClicked(int)),this,SLOT(setDrawMode1D(int)));
    tool_bar_->addSeparator();


    link_box_ = new QComboBox(tool_bar_1d_);
    link_box_->setToolTip("Linking spectra");
    link_box_->setWhatsThis("Linking spectra<BR><BR>Use this combobox to link two 1D spectra."
    												"Linked spectra zoom in/out together.");
    tool_bar_1d_->addWidget(link_box_);
    connect(link_box_,SIGNAL(activated(int)),this,SLOT(linkActiveTo(int)));

    //--2D toolbar--
    tool_bar_2d_ = addToolBar("2D tool bar");

    dm_precursors_2d_ = tool_bar_2d_->addAction(QPixmap(precursors),"Show MS/MS precursors");
    dm_precursors_2d_->setCheckable(true);
    dm_precursors_2d_->setWhatsThis("2D peak draw mode: Precursors<BR><BR>MS/MS precursor peaks are marked");
    connect(dm_precursors_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    projections_2d_ = tool_bar_2d_->addAction(QPixmap(projections), "Show Projections" ,this, SLOT(showProjections()));
    projections_2d_->setWhatsThis("Projections: Shows projections of peak data along RT and MZ axis.");


    dm_hull_2d_ = tool_bar_2d_->addAction(QPixmap(convexhulls),"Show feature convex hulls");
    dm_hull_2d_->setCheckable(true);
    dm_hull_2d_->setWhatsThis("2D feature draw mode: Convex hull<BR><BR>The convex hull of the feature is displayed");
    connect(dm_hull_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    dm_numbers_2d_ = tool_bar_2d_->addAction(QPixmap(numbers),"Show feature numbers");
    dm_numbers_2d_->setCheckable(true);
    dm_numbers_2d_->setWhatsThis("2D feature draw mode: Numbers<BR><BR>The feature number is displayed next to the feature");
    connect(dm_numbers_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    //layer wndow
    QDockWidget* layer_bar = new QDockWidget("Layers", this);
    addDockWidget(Qt::RightDockWidgetArea, layer_bar);
    layer_manager_ = new QListWidget(layer_bar);
    layer_manager_->setWhatsThis("Layer bar<BR><BR>Here the availabe layers are shown. You can select, hide"
    								              " and remove the layers using this bar.");

    layer_bar->setWidget(layer_manager_);
    layer_manager_->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(layer_manager_,SIGNAL(currentRowChanged(int)),this,SLOT(layerSelectionChange(int)));
		connect(layer_manager_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(layerContextMenu(const QPoint&)));
		connect(layer_manager_,SIGNAL(itemChanged(QListWidgetItem*)),this,SLOT(layerVisibilityChange(QListWidgetItem*)));
    windows->addAction("&Show layer window",layer_bar,SLOT(show()));

		//log window
		QDockWidget* log_bar = new QDockWidget("Log", this);
		addDockWidget(Qt::BottomDockWidgetArea, log_bar);
		log_ = new QTextEdit(log_bar);
		log_->setReadOnly(true);
		log_bar->setWidget(log_);
		log_bar->hide();
    windows->addAction("&Show log window",log_bar,SLOT(show()));

	//################## DEFAULTS #################
    //general
    defaults_.setValue("preferences:default_map_view", "2d","Default visualization mode for maps.");
    defaults_.setValue("preferences:default_path", ".","Default path for loading and storing files.");
    defaults_.setValue("preferences:tmp_file_path", "/tmp/","Path where temporary files can be created.");
    defaults_.setValue("preferences:number_of_recent_files", 15,"Number of recent files in the main menu.");
    defaults_.setValue("preferences:legend", "show", "Legend visibility ('show' or 'hide')");
    defaults_.setValue("preferences:intensity_cutoff", "none","Low intensity cutoff for maps.");
    //db
    defaults_.setValue("preferences:db:host", "localhost", "Database server host name.");
    defaults_.setValue("preferences:db:login", "NoName", "Database login.");
    defaults_.setValue("preferences:db:name", "OpenMS", "Database name.");
    defaults_.setValue("preferences:db:port", "3306", "Database server port.");
    defaults_.setSectionDescription("preferences:db","Database settings.");
    //1d
    Spectrum1DCanvas* def1 = new Spectrum1DCanvas(Param(),0);
    defaults_.insert("preferences:1d:",def1->getDefaults());
    delete def1;
    defaults_.setSectionDescription("preferences:1d","Settings for single spectrum view.");
    //2d
    Spectrum2DCanvas* def2 = new Spectrum2DCanvas(Param(),0);
    defaults_.insert("preferences:2d:",def2->getDefaults());
    defaults_.setSectionDescription("preferences:2d","Settings for 2D map view.");
    delete def2;
    //3d
    Spectrum3DCanvas* def3 = new Spectrum3DCanvas(Param(),0);
    defaults_.insert("preferences:3d:",def3->getDefaults());
    delete def3;
    defaults_.setSectionDescription("preferences:3d","Settings for 3D map view.");
		
		defaults_.setValue("preferences:version","none","OpenMS version, used to check if the TOPPView.ini is up-to-date");
		
		subsections_.push_back("preferences:RecentFiles");
		
  	defaultsToParam_();
  	
  	//load param file
    loadPreferences();    
  }

  TOPPViewBase::~TOPPViewBase()
  {
  	//cout << "DEST TOPPViewBase" << endl;
  	savePreferences();
  }

  void TOPPViewBase::closeEvent(QCloseEvent* /*event*/)
  {
  	ws_->closeAllWindows();
  }

  void TOPPViewBase::addDBSpectrum(UInt db_id, bool as_new_window, bool maps_as_2d, bool maximize, OpenDialog::Mower use_mower)
  {
    //DBConnection for all DB queries
    DBConnection con;
    con.connect(param_.getValue("preferences:db:name"), param_.getValue("preferences:db:login"),param_.getValue("DBPassword"),param_.getValue("preferences:db:host"),(Int)param_.getValue("preferences:db:port"));

    //DB adapter
    DBAdapter dba(con);

    String db_id_string(db_id);
    QSqlQuery result;
    con.executeQuery("SELECT count(id) from DATA_Spectrum where fid_MSExperiment='"+db_id_string+"' and MSLevel='1'",result);

    //tab caption
    String caption = "DB ("+db_id_string+")";

    SpectrumWidget* w;

    SpectrumCanvas::ExperimentType* exp = 0;

    if (activeWindow_()==0)
    {
      as_new_window = true;
    }

    //open in new window
    if (as_new_window)
    {
      //create 1D View
      if (result.value(0).toInt()==1)
      {
      	//cout << "NEW 1D" << endl;
        // create 1D window
        w = new Spectrum1DWidget(param_.copy("preferences:1d:",true), ws_);

        //determine Spectrum id
        con.executeQuery("SELECT id from DATA_Spectrum where fid_MSExperiment='"+db_id_string+"' and MSLevel='1'",result);
        UID spectrum_id = result.value(0).toInt();

        //load data
        exp = &(w->canvas()->addEmptyPeakLayer());
        exp->resize(1);
        dba.loadSpectrum(spectrum_id,(*exp)[0]);
      }
      //create 2D/3D view
      else
      {
        //create 2D view
        if (maps_as_2d)
        {
          //cout << "NEW 2D" << endl;
          //create 2D window
          w = new Spectrum2DWidget(param_.copy("preferences:2d:",true), ws_);

          //load spectrum
          exp = &(w->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, (*exp));
        }
        //create 3D view
        else
        {
          //cout << "NEW 3D" << endl;
        	// create 3D window
          w = new Spectrum3DWidget(param_.copy("preferences:3d:",true), ws_);

          //load data
          exp = &(w->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, *exp);
        }
      }
    }
    //open in active window
    else
    {
      //create 1D View
      if (result.value(0).toInt()==1)
      {
      	//cout << "ACTIVE 1D" << endl;
        w = active1DWindow_();
        //wrong active window type
        if (w==0)
        {
          QMessageBox::critical(this,"Wrong file type",("You cannot open 1D data ("+db_id_string+") in a 2D window!<BR>Please open the file in new tab.").c_str());
          return;
        }
        else //open it
        {
          //determine Spectrum id
          con.executeQuery("SELECT id from DATA_Spectrum where fid_MSExperiment='"+db_id_string+"' and MSLevel='1'", result);
          UID spectrum_id = result.value(0).toInt();

          //load data
	        exp = &(w->canvas()->addEmptyPeakLayer());
	        exp->resize(1);
	        dba.loadSpectrum(spectrum_id,(*exp)[0]);
        }
      }
      //create 2D/3D view
      else
      {
        Spectrum2DWidget* w2 = active2DWindow_();
        Spectrum3DWidget* w3 = active3DWindow_();
        //wrong active window type
        if (w2==0 && w3==0)
        {
          QMessageBox::critical(this,"Wrong file type",("You cannot open 1D data ("+db_id_string+") in a 2D/3D window!<BR>Please open the file in new tab.").c_str());
          return;
        }
        //create 2D view
        if (w2!=0)
        {
        	//cout << "ACTIVE 2D" << endl;
          w = w2;

          //load data
          exp = &(w->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, (*exp));
        }
        //create 3D view
        else
        {
        	//cout << "ACTIVE 3D" << endl;
          w = w3;

          //load data
          exp = &(w->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, (*exp));
        }
      }
		}
		
    //noise estimator
    float cutoff = 0;
    if(use_mower!=OpenDialog::NO_MOWER && exp->size()>1)
    {
      cutoff = estimateNoise_(*exp);
    }
    w->canvas()->finishAdding(cutoff);
		w->canvas()->setCurrentLayerName(caption);
    //use_mower

    //do for all windows
    if (as_new_window)
    {
      showAsWindow_(w,caption);
    }

    //do for all (in active and in new window, 1D/2D/3D)
    if(maximize)
    {
      w->showMaximized();
    }

    //do for all windows
    updateLayerbar();
  }

  float TOPPViewBase::estimateNoise_(const SpectrumCanvas::ExperimentType& exp)
  {
  	//test if no scans with MS-level 1 exist => prevent deadlock
  	bool ms1_present = false;
  	for (UInt i = 0; i < exp.size(); ++i)
  	{
  		if (exp[i].getMSLevel()==1)
  		{
  			ms1_present = true;
  			break;
  		}
  	}
  	if (!ms1_present)
  	{
  		return 0.0;
  	}
  	
    float noise = 0.0;
    UInt count = 0;
    srand(time(0));
    //cout << "size: " << exp.size() << endl;
    while (count<10)
    {
      UInt scan = (UInt)( (double)rand() / ((double)(RAND_MAX)+1.0f) * (double)(exp.size()-1) );

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
        //cout << "scan: "<< scan <<" Groesse: " << tmp.size() << " Index: " << (UInt)ceil((float)(tmp.size()-1)/1.25f) << " Wert: "<< tmp[(UInt)ceil((float)(tmp.size()-1)/1.25f)] << endl;
        noise += tmp[(UInt)ceil((float)(tmp.size()-1)/1.25f)];

        ++count;
      }
    }
    noise = noise / 10.0f;
    return noise;
  }

  void TOPPViewBase::preferencesDialog()
  {
		Internal::TOPPViewPrefDialog dlg(this);
		
		//get pointers
		QLineEdit* default_path = dlg.findChild<QLineEdit*>("default_path");
		QLineEdit* temp_path = dlg.findChild<QLineEdit*>("temp_path");
		QSpinBox* recent_files = dlg.findChild<QSpinBox*>("recent_files");
		QComboBox* map_default = dlg.findChild<QComboBox*>("map_default");
		QComboBox* map_cutoff = dlg.findChild<QComboBox*>("map_cutoff");
		
		QLineEdit* db_host = dlg.findChild<QLineEdit*>("db_host");
		QSpinBox* db_port = dlg.findChild<QSpinBox*>("db_port");
		QLineEdit* db_name = dlg.findChild<QLineEdit*>("db_name");
		QLineEdit* db_login = dlg.findChild<QLineEdit*>("db_login");
		
		ColorSelector* color_1D = dlg.findChild<ColorSelector*>("color_1D");
		ColorSelector* selected_1D = dlg.findChild<ColorSelector*>("selected_1D");
		ColorSelector* icon_1D = dlg.findChild<ColorSelector*>("icon_1D");

		MultiGradientSelector* peak_2D = dlg.findChild<MultiGradientSelector*>("peak_2D");
		QComboBox* mapping_2D = dlg.findChild<QComboBox*>("mapping_2D");

		MultiGradientSelector* peak_3D = dlg.findChild<MultiGradientSelector*>("peak_3D");
		QComboBox* shade_3D = dlg.findChild<QComboBox*>("shade_3D");
		QSpinBox* line_width_3D  = dlg.findChild<QSpinBox*>("line_width_3D");
		QComboBox* reduction_3D = dlg.findChild<QComboBox*>("reduction_3D");
		QSpinBox* reduction_peaks_3D  = dlg.findChild<QSpinBox*>("reduction_peaks_3D");
		
		//set General values
		default_path->setText(param_.getValue("preferences:default_path").toQString());
		temp_path->setText(param_.getValue("preferences:tmp_file_path").toQString());
		recent_files->setValue((Int)param_.getValue("preferences:number_of_recent_files"));
		map_default->setCurrentIndex(map_default->findText(param_.getValue("preferences:default_map_view").toQString()));
		map_cutoff->setCurrentIndex(map_cutoff->findText(param_.getValue("preferences:intensity_cutoff").toQString()));		

		db_host->setText(param_.getValue("preferences:db:host").toQString());
		db_port->setValue((Int)param_.getValue("preferences:db:port"));
		db_name->setText(param_.getValue("preferences:db:name").toQString());
		db_login->setText(param_.getValue("preferences:db:login").toQString());
		
		color_1D->setColor(QColor(param_.getValue("preferences:1d:peak_color").toQString()));
		selected_1D->setColor(QColor(param_.getValue("preferences:1d:highlighted_peak_color").toQString()));
		icon_1D->setColor(QColor(param_.getValue("preferences:1d:icon_color").toQString()));

		peak_2D->gradient().fromString(param_.getValue("preferences:2d:dot:gradient"));
		mapping_2D->setCurrentIndex(mapping_2D->findText(param_.getValue("preferences:2d:mapping_of_mz_to").toQString()));

		peak_3D->gradient().fromString(param_.getValue("preferences:3d:dot:gradient"));
		shade_3D->setCurrentIndex((Int)param_.getValue("preferences:3d:dot:shade_mode"));
		line_width_3D->setValue((Int)param_.getValue("preferences:3d:dot:line_width"));
		reduction_3D->setCurrentIndex(reduction_3D->findText(param_.getValue("preferences:3d:reduction_mode").toQString()));	
		reduction_peaks_3D->setValue((Int)param_.getValue("preferences:3d:displayed_peaks"));
		
		//execute dialog
		if (dlg.exec())
		{
			//load data to param
			param_.setValue("preferences:default_path", default_path->text().toAscii().data());
			param_.setValue("preferences:tmp_file_path", temp_path->text().toAscii().data());
			param_.setValue("preferences:number_of_recent_files", recent_files->value());
			param_.setValue("preferences:default_map_view", map_default->currentText().toAscii().data());
			param_.setValue("preferences:intensity_cutoff", map_cutoff->currentText().toAscii().data());

			param_.setValue("preferences:db:host",db_host->text().toAscii().data());
			param_.setValue("preferences:db:port",db_port->value());
			param_.setValue("preferences:db:name",db_name->text().toAscii().data());
			param_.setValue("preferences:db:login",db_login->text().toAscii().data());
			param_.remove("DBPassword");

			param_.setValue("preferences:1d:peak_color",color_1D->getColor().name().toAscii().data());
			param_.setValue("preferences:1d:highlighted_peak_color",selected_1D->getColor().name().toAscii().data());
			param_.setValue("preferences:1d:icon_color",icon_1D->getColor().name().toAscii().data());

			param_.setValue("preferences:2d:dot:gradient",peak_2D->gradient().toString());
			param_.setValue("preferences:2d:mapping_of_mz_to",mapping_2D->currentText().toAscii().data());

			param_.setValue("preferences:3d:dot:gradient",peak_3D->gradient().toString());
			param_.setValue("preferences:3d:dot:shade_mode", shade_3D->currentIndex());
			param_.setValue("preferences:3d:dot:line_width",line_width_3D->value());
			param_.setValue("preferences:3d:reduction_mode", reduction_3D->currentText().toAscii().data());
			param_.setValue("preferences:3d:displayed_peaks",	reduction_peaks_3D->value());

			savePreferences();
		}
  }

  void TOPPViewBase::addSpectrum(const String& filename,bool as_new_window, bool maps_as_2d, bool maximize, OpenDialog::Mower use_mower, FileHandler::Type force_type)
  {
    if (!File::exists(filename))
    {
      QMessageBox::critical(this,"Open file error",("The file '"+filename+"' does not exist!").c_str());
      return;
    }

    //extract the filename without path
    String caption = File::basename(filename);

    //update recent files list
    addRecentFile_(filename);

    //windowpointer
    SpectrumWidget* w=0;
    
    if (activeWindow_()==0)
    {
      as_new_window = true;
    }

		//determine file type if not forced
		FileHandler fh;
		if (force_type==FileHandler::UNKNOWN)
		{
			if (force_type==FileHandler::UNKNOWN)
			{
				force_type = fh.getTypeByFileName(filename);
			}

			if (force_type==FileHandler::UNKNOWN)
			{
				force_type = fh.getTypeByContent(filename);
			}
		}

		if (force_type==FileHandler::UNKNOWN)
		{
      QMessageBox::critical(this,"Open file error",("Could not determine file type of '"+filename+"'!").c_str());
      return;
		}

#ifdef DEBUG_TOPP
		cout << "TOPPViewBase::addSpectrum():";
		cout << " - File name: " << filename << endl;
		cout << " - File type: " << fh.typeToName(force_type) << endl;
		cout << " - New Window: " << as_new_window << endl;
		cout << " - Map as 2D: " << maps_as_2d << endl;		
#endif

    if (as_new_window)
    {
      if (force_type==FileHandler::DTA)
      {
        w = new Spectrum1DWidget(param_.copy("preferences:1d:",true), ws_);
      }
      else if (maps_as_2d || force_type==FileHandler::FEATURE || force_type==FileHandler::FEATURE_PAIRS) //2d or features
      {
        w = new Spectrum2DWidget(param_.copy("preferences:2d:",true), ws_);
      }
      else //3d
      {
        w = new Spectrum3DWidget(param_.copy("preferences:3d:",true), ws_);
      }
    }
    else //!as_new_window
    {
      if (active1DWindow_()!=0) //1D window
      {
        w = active1DWindow_();
      }
      else if (active2DWindow_()!=0) //2d window
      {
        if (force_type==FileHandler::DTA)
        {
          QMessageBox::critical(this,"Wrong file type",("You cannot open 1D data ("+filename+") in a 2D window!").c_str());
          return;
        }
        w = active2DWindow_();
      }
      else if (active3DWindow_()!=0)//3d window
      {
        w = active3DWindow_();
      }
    }

    //try to read the data from file
    if (force_type==FileHandler::FEATURE) //features
    {
      FeatureMap<> map;
      try
      {
        FeatureXMLFile().load(filename,map);
      }
      catch(Exception::Base& e)
      {
        QMessageBox::critical(this,"Error",(String("Error while reading feature file: ")+e.what()).c_str());
        return;
      }
      w->canvas()->addLayer(map,false);
      w->canvas()->setCurrentLayerName(caption);
    }
    else if (force_type==FileHandler::FEATURE_PAIRS) //feature pairs
    {
    	//load pairs
      std::vector< ElementPair < Feature > >  pairs;
      try
      {
        FeaturePairsXMLFile().load(filename,pairs);
      }
      catch(Exception::Base& e)
      {
        QMessageBox::critical(this,"Error",(String("Error while reading feature pairs file: ")+e.what()).c_str());
        return;
      }
      
      //convert to features
      FeatureMap<> map;
      FeaturePairsXMLFile::pairsToFeatures(pairs,map);
      w->canvas()->addLayer(map,true);
      w->canvas()->setCurrentLayerName(caption);
    }
    else
    {
      //try to read the data from file (raw/peak data)
      SpectrumCanvas::ExperimentType* exp = &(w->canvas()->addEmptyPeakLayer());
      try
      {
        fh.loadExperiment(filename,*exp, force_type,ProgressLogger::GUI);
      }
      catch(Exception::Base& e)
      {
        QMessageBox::critical(this,"Error",(String("Error while reading data file: ")+e.what()).c_str());
        return;
      }

      //check if only one scan is in a 2d file
      if (as_new_window && active1DWindow_()==0 && exp->size()==1)
      {
        delete(w);
        w = new Spectrum1DWidget(param_.copy("preferences:1d:",true), ws_);
        exp = &(w->canvas()->addEmptyPeakLayer());
        FileHandler().loadExperiment(filename,*exp, force_type);
      }

      //do for all (in active and in new window, 1D/2D/3D)
      float cutoff=0;

      //calculate noise
      if(use_mower!=OpenDialog::NO_MOWER && exp->size()>1)
      {
        cutoff = estimateNoise_(*exp);
      }
      w->canvas()->finishAdding(cutoff);
      w->canvas()->setCurrentLayerName(caption);
    }
    
  	updateLayerbar();

    if (as_new_window)
    {
      showAsWindow_(w,caption);
    }
    if(maximize)
    {
      w->showMaximized();
    }
  }

  void TOPPViewBase::addRecentFile_(const String& filename)
  {
    String tmp = filename;

    //add prefix to relative paths
    if (!tmp.hasPrefix("/"))
    {
      tmp = QDir::currentPath().toAscii().data()+string("/")+ tmp;
    }
	
		// remove the new file if already in the recent list and prepend it
		recent_files_.removeAll(tmp.c_str());
		recent_files_.prepend(tmp.c_str());
		
		//remove those files exceeding the defined number
		UInt number_of_recent_files = UInt(param_.getValue("preferences:number_of_recent_files"));
		while ((UInt)recent_files_.size() > number_of_recent_files)
		{
			recent_files_.removeLast();
		}
		
    updateRecentMenu_();
  }

  void TOPPViewBase::updateRecentMenu_()
  {
    //get/correct number of recent files
		UInt number_of_recent_files = UInt(param_.getValue("preferences:number_of_recent_files"));
		if (number_of_recent_files>20)
		{
			number_of_recent_files = 20;
			param_.setValue("preferences:number_of_recent_files",20);
		}
		
		for (UInt i = 0; i < 20; ++i)
		{
			if (i < (UInt)(recent_files_.size()))
			{
				recent_actions_[i]->setText(recent_files_[i]);
				recent_actions_[i]->setVisible(true);
			}
			else
			{
				recent_actions_[i]->setVisible(false);
			}
		}
  }

  void TOPPViewBase::maximizeActiveSpectrum()
  {
    if (ws_->activeWindow())
    {
      ws_->activeWindow()->showMaximized();
    }
  }

  void TOPPViewBase::addTab_(SpectrumWidget* w, const String& tabCaption)
  {
  	//static window counter
  	static int window_counter = 4711;
  	w->window_id = ++window_counter;

		//add tab and assign window id as data  	
    int tab_index = tab_bar_->addTab(tabCaption.c_str());
    tab_bar_->setTabData(tab_index, window_counter);
    //cout << "Added tab: '" << tabCaption << "' => " << window_counter << endl;
        
    //connect slots and sigals for removing the spectrum from the bar, when it is closed
    connect(w,SIGNAL(aboutToBeDestroyed(int)),this,SLOT(removeTab(int)));
    tab_bar_->setCurrentIndex(tab_index);
  }

  SpectrumWidget* TOPPViewBase::window_(int id) const
  {
  	//cout << "Looking for tab with id: " << id << endl;
  	QList<QWidget*> windows = ws_->windowList();
		for(int i=0; i< windows.size(); ++i)
		{
			SpectrumWidget* window = dynamic_cast<SpectrumWidget*>(windows.at(i));
			//cout << "  Tab " << i << ": " << window->window_id << endl;
			if (window->window_id == id)
			{
				return window;
			}
		}
		return 0;
  }

  void TOPPViewBase::removeTab(int id)
  {
  	for (int i=0; i < tab_bar_->count(); ++i)
  	{ 
  		if( tab_bar_->tabData(i).toInt() == id)
  		{
  			tab_bar_->removeTab(i);
  			return;
  		}
  	}
  }

  void TOPPViewBase::closeByTab(int index)
  {
  	SpectrumWidget* window = window_(tab_bar_->tabData(index).toInt());
  	if (window)
  	{
  		window->close();
  	}
  }
 
  void TOPPViewBase::focusByTab(int index)
  {
  	SpectrumWidget* window = window_(tab_bar_->tabData(index).toInt());
  	if (window)
  	{
  		window->setFocus();
  	}
  }

  void TOPPViewBase::saveImage()
  {
    //check if there is a active window
    SpectrumWidget* window = activeWindow_();
    if (window!=0)
    {
      SaveImageDialog* dialog = new SaveImageDialog(this);
      dialog->setSize(1024,768);
      if (dialog->exec())
      {
        QString format=dialog->getFormat();
        QString file_name = QFileDialog::getSaveFileName(this, "Save file", window->windowTitle().section('.',0,0),
                            format+" Images(*."+format+" *."+format.toLower()+")");
        if (!file_name.isEmpty())
        {
          //append missing file extension
          if (!file_name.toUpper().endsWith(format))
          {
            file_name.append("."+format.toLower());
          }
					if (File::writable(file_name))
					{
          	QImage image = window->getImage(dialog->getXSize(),dialog->getYSize());
          	image.save(file_name,format.toAscii().data(),100);
        	}
        	else
        	{
        	  QMessageBox::critical(this,"Error writing file!",(String("Cannot write to '")+file_name	+"'!").c_str());
        	}
        }
      }
    }
  }

  void TOPPViewBase::print()
  {
#ifndef QT_NO_PRINTER
    SpectrumWidget* window = activeWindow_();
    if (window!=0)
    {
      QPrinter* printer = new QPrinter(QPrinter::HighResolution);
      printer->setResolution(300);
      
      QPrintDialog dialog(printer, this);
 			if (dialog.exec())
      {
        QPainter p;
        if (!p.begin(printer)) return;
        unsigned int dpix = printer->logicalDpiX();
        unsigned int dpiy = printer->logicalDpiY();
        QRect body(dpix,dpiy,printer->width()-2*dpix,printer->height()-2*dpiy);			// one inch margin
        QImage image = window->getImage(body.width(),body.height());
        p.drawImage(body,image);
        QString titel = QString("%1\n%2").arg(window->windowTitle().section('/',-1)).arg(QDate::currentDate().toString());
        //	p.drawText(dpix,0,body.width(),dpiy, Qt::AlignCenter, titel);
        p.drawText(dpix,body.height()+dpiy, body.width(), dpiy, Qt::AlignCenter, titel);
      }
      delete(printer);
    }
#endif

  }

  void TOPPViewBase::closeFile()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      ws_->activeWindow()->close();
    }
  }

  void TOPPViewBase::saveLayer()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
		
      //warn if hidden layer => wrong layer selected...
    	if (!layer.visible)
    	{
    		QMessageBox::warning(this,"Warning","The current layer is not visible!");
    	}
    	//Visible area
    	DoubleReal min_rt;
    	DoubleReal max_rt;
    	DoubleReal min_mz;
    	DoubleReal max_mz;
    	
			if (dynamic_cast<Spectrum3DCanvas*>(activeWindow_()->canvas())!=0) //3D
			{
		    	min_rt = activeWindow_()->canvas()->getVisibleArea().min()[0];
		    	max_rt = activeWindow_()->canvas()->getVisibleArea().max()[0];
		    	min_mz = activeWindow_()->canvas()->getVisibleArea().min()[1];
		    	max_mz = activeWindow_()->canvas()->getVisibleArea().max()[1];
			}
			else //1D or 2D
			{
		    	min_rt = activeWindow_()->canvas()->getVisibleArea().min()[1];
		    	max_rt = activeWindow_()->canvas()->getVisibleArea().max()[1];
		    	min_mz = activeWindow_()->canvas()->getVisibleArea().min()[0];
		    	max_mz = activeWindow_()->canvas()->getVisibleArea().max()[0];			
			}
    	
    	//cout << "RT: " << min_rt << "-"  << max_rt << " -- mz: " << min_mz << "-" << max_mz << endl; 
    	
    	if (layer.type==LayerData::DT_PEAK)
    	{
    		//Extract selected visible data to out
    		LayerData::ExperimentType out;
    		out.ExperimentalSettings::operator=(layer.peaks);
    		LayerData::ExperimentType::ConstIterator begin;
    		LayerData::ExperimentType::ConstIterator end; 
    		if (layer.peaks.size()==1)
    		{
	    		begin = layer.peaks.begin();
	    		end = layer.peaks.end();
    		}
    		else
    		{
	    		begin = layer.peaks.RTBegin(min_rt);
	    		end = layer.peaks.RTEnd(max_rt);
    		}
    		out.resize(end-begin);
				
				UInt i = 0;
    		for (LayerData::ExperimentType::ConstIterator it=begin; it!=end; ++it)
    		{
  				out[i].SpectrumSettings::operator=(*it);
  				out[i].setRT(it->getRT());
  				out[i].setMSLevel(it->getMSLevel());
  				out[i].setPrecursorPeak(it->getPrecursorPeak());
  				for (LayerData::ExperimentType::SpectrumType::ConstIterator it2 = it->MZBegin(min_mz); it2!= it->MZEnd(max_mz); ++it2)
  				{
  					if ( it2->getIntensity() >= layer.min_int && it2->getIntensity() <= layer.max_int)
  					{
  						out[i].push_back(*it2);
  					}
  				}
  				++i;
    		}
    		//no extracted data
    		if (out.size()==0)
    		{
    		  QMessageBox::critical(this,"Error","The displayed region of the current layer is empty!");
    		  return;
    		}
    		//one scan => DTA
    		else if (out.size()==1)
    		{
		      QString file_name = QFileDialog::getSaveFileName(this, "Save file", param_.getValue("preferences:default_path").toString().c_str(),
					                    "DTA files (*.dta)");
					if (!file_name.isEmpty())
					{
					  DTAFile().store(file_name.toAscii().data(),out[0]);
					}
    		}
    		//more than one scan => MzData
    		else
    		{
		      QString file_name = QFileDialog::getSaveFileName(this, "Save file",  param_.getValue("preferences:default_path").toString().c_str(),
					                    "MzData files (*.mzData)");
					if (!file_name.isEmpty())
					{
					  MzDataFile f;
					  f.setLogType(ProgressLogger::GUI);
					  f.store(file_name.toAscii().data(),out);
					}
    		}
			}
			else //feature data
			{
    		//Extract selected visible data to out
    		LayerData::FeatureMapType out;
    		out.ExperimentalSettings::operator=(layer.features);
    		for (LayerData::FeatureMapType::ConstIterator it=layer.features.begin(); it!=layer.features.end(); ++it)
    		{
					if ( it->getIntensity() >= layer.min_int && 
							 it->getIntensity() <= layer.max_int &&
							 it->getRT() >= min_rt &&
							 it->getRT() <= max_rt &&
							 it->getMZ() >= min_mz &&
							 it->getMZ() <= max_mz
						 )
					{
						out.push_back(*it);
					}
  			}
    		//no extracted data
    		if (out.size()==0)
    		{
    		  QMessageBox::critical(this,"Error","The displayed region of the current layer is empty!");
    		  return;
    		}
    		else
    		{
		      QString file_name = QFileDialog::getSaveFileName(this, "Save file",  param_.getValue("preferences:default_path").toString().c_str(),
					                    "features files (*.featureXML)" );
					if (!file_name.isEmpty())
					{
					  FeatureXMLFile().store(file_name.toAscii().data(),out);
					}
    		}
			}
    }
  }

  void TOPPViewBase::editMetadata()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
      //warn if hidden layer => wrong layer selected...
    	if (!layer.visible)
    	{
    		QMessageBox::warning(this,"Warning","The current layer is not visible!");
    	}
			MSMetaDataExplorer dlg(true, this);
      dlg.setWindowTitle("Edit meta data");
			if (layer.type==LayerData::DT_PEAK) //peak data
    	{
    		dlg.visualize(const_cast<LayerData&>(layer).peaks);
				
    	}
    	else //feature data
    	{
    		dlg.visualize(const_cast<LayerData&>(layer).features);
    	}
      dlg.exec();
    }
  }

  void TOPPViewBase::layerPreferencesDialog()
  {
    if (ws_->activeWindow())
    {
			activeWindow_()->canvas()->showCurrentLayerPreferences();
    }
  }

  void TOPPViewBase::layerIntensityDistribution()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      activeWindow_()->showIntensityDistribution();
    }  	
  }

  void TOPPViewBase::changeAxisVisibility()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      activeWindow_()->showLegend(!activeWindow_()->isLegendShown());
    }  	
  }

  void TOPPViewBase::linkActiveTo(int index)
  {
  	Spectrum1DWidget* active = active1DWindow_();
  	if (active==0) return;
  	
  	//cout << "linkActiveTo() active: " << active->window_id << endl;
  	
  	//remove link if present
  	if (link_map_.find(active->window_id)!=link_map_.end())
  	{
  		SpectrumWidget* linked_to = window_(link_map_[active->window_id]);
  		if (linked_to != 0)
  		{
  			//cout << "  disconnect signals to: " << linked_to->window_id << endl;
	  		//disconnect outgoing signals
	  		disconnect(active->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),linked_to->canvas(),SLOT(setVisibleArea(DRange<2>)));
  			//disconnect incoming signals
  			disconnect(linked_to->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),active->canvas(),SLOT(setVisibleArea(DRange<2>)));
  			//erase entries in id map
  			link_map_.erase(active->window_id);
  			link_map_.erase(linked_to->window_id);
  		}
  		else
  		{
  			cout << "linkActiveTo() - disconnect: Error, could not find window with id '" << active->window_id << "'!" << endl;
  		}
  	}
  	  
    if (link_box_->itemText(index)=="<unlinked>") return;
    
    //link
    SpectrumWidget* link_to = window_(link_box_->itemData(index).toInt());
		if (link_to!=0)
		{
			//cout << "  connecting signals to: " << link_to->window_id << endl;
		  connect(active->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),link_to->canvas(),SLOT(setVisibleArea(DRange<2>)));
      connect(link_to->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),active->canvas(),SLOT(setVisibleArea(DRange<2>)));
      //add entries to id map
      link_map_[link_to->window_id]=active->window_id;
      link_map_[active->window_id]=link_to->window_id;
		}
		else
		{
			cout << "linkActiveTo() - connect: Error, could not find window with id '" << link_box_->itemData(index).toInt() << "'!" << endl;
		}    
  }

  void TOPPViewBase::showStatusMessage(string msg, OpenMS::UInt time)
  {
    if (time==0)
    {
      message_label_->setText(msg.c_str());
      statusBar()->update();
    }
    else
    {
      statusBar()->showMessage(msg.c_str(), time);
    }
  }

  void TOPPViewBase::showCursorStatus(double mz, double intensity, double rt)
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
    statusBar()->update();
  }

  void TOPPViewBase::changeGridLines()
  {
    SpectrumWidget* window = activeWindow_();
    if (window!=0)
    {
      window->canvas()->showGridLines(!window->canvas()->gridLinesShown());
    }
  }

  void TOPPViewBase::resetZoom()
  {
    SpectrumWidget* window = activeWindow_();
    if (window!=0)
    {
      window->canvas()->resetZoom();
    }
  }

  void TOPPViewBase::setActionMode(int index)
  {
    SpectrumWidget* w = activeWindow_();
    if (w)
    {
    	w->setActionMode((OpenMS::SpectrumCanvas::ActionModes)index);
  	}
  }

  void TOPPViewBase::setIntensityMode(int index)
  {
    SpectrumWidget* w = activeWindow_();
    if (w)
    {
    	w->setIntensityMode((OpenMS::SpectrumCanvas::IntensityModes)index);
  	}
  }

  void TOPPViewBase::setDrawMode1D(int index)
  {
    Spectrum1DWidget* w = active1DWindow_();
    if (w)
    {
    	w->canvas()->setDrawMode((OpenMS::Spectrum1DCanvas::DrawModes)index);
  	}
  }

  void TOPPViewBase::changeLayerFlag(bool on)
  {
		QAction* action = qobject_cast<QAction *>(sender());
    if (Spectrum2DWidget* win = active2DWindow_())
    {
    	//peaks
			if (action == dm_precursors_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::P_PRECURSORS,on);
			}
			//features
			else if (action == dm_hull_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::F_HULLS,on);
			}
			else if (action == dm_numbers_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::F_NUMBERS,on);
			}
		}
  }

  void TOPPViewBase::updateToolbar()
  {
    SpectrumWidget* w = activeWindow_();

    if (w)
    {
      //set action mode
     	action_group_->button(w->getActionMode())->setChecked(true);

      //set intensity mode
     	intensity_group_->button(w->canvas()->getIntensityMode())->setChecked(true);
    }

    //1D
    Spectrum1DWidget* w1 = active1DWindow_();
    if (w1)
    {
      //draw mode
      draw_group_1d_->button(w1->canvas()->getDrawMode())->setChecked(true);

      //update link selector
      int item_index = -1;
      link_box_->clear();
      link_box_->insertItem(++item_index,"<unlinked>",0);
      QWidgetList windows = ws_->windowList();
      for ( int i = 0; i < windows.count(); ++i )
      {
        Spectrum1DWidget* window = dynamic_cast<Spectrum1DWidget*>(windows.at(i));
        if (window !=0 && window!=w)
        {
          link_box_->insertItem(++item_index,File::basename(window->windowTitle().toAscii().data()).c_str(),window->window_id);
        	if (link_map_.find(w1->window_id)!=link_map_.end() && link_map_[w1->window_id] == window->window_id)
          {
            link_box_->setCurrentIndex(item_index);
          }
        }
      }

      //show/hide toolbars and buttons
      tool_bar_1d_->show();
      tool_bar_2d_->hide();
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_SELECT)->setEnabled(true);
    }

    //2d
    Spectrum2DWidget* w2 = active2DWindow_();
    if (w2)
    {
      //peak draw modes
      if (w2->canvas()->getCurrentLayer().type == LayerData::DT_PEAK)
      {
      	dm_precursors_2d_->setVisible(true);
      	projections_2d_->setVisible(true);
      	dm_hull_2d_->setVisible(false);
      	dm_numbers_2d_->setVisible(false);
				dm_precursors_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::P_PRECURSORS));
			}
			//feature draw modes
			else
			{
      	dm_precursors_2d_->setVisible(false);
      	projections_2d_->setVisible(false);
      	dm_hull_2d_->setVisible(true);
      	dm_numbers_2d_->setVisible(true);
      	dm_hull_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::F_HULLS));
      	dm_numbers_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::F_NUMBERS));
			}
      //show/hide toolbars and buttons
      tool_bar_1d_->hide();
      tool_bar_2d_->show();
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_SELECT)->setEnabled(true);
    }

    //1D
    Spectrum3DWidget* w3 = active3DWindow_();
    if (w3)
    {
      //show/hide toolbars and buttons
      tool_bar_1d_->hide();
      tool_bar_2d_->hide();
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_SELECT)->setEnabled(false);
    }
  }

  void TOPPViewBase::updateLayerbar()
  {
		layer_manager_->clear();

    SpectrumCanvas* cc = activeCanvas_();
    if (cc == 0)
    {
      return;
    }
		QListWidgetItem* item = 0;
    for (UInt i = 0; i<cc->getLayerCount(); ++i)
    {
    	//add item
    	item = new QListWidgetItem( layer_manager_ );
			item->setText(cc->getLayer(i).name.c_str());
    	if (cc->getLayer(i).visible)
    	{
    		item->setCheckState(Qt::Unchecked);
    	}
    	else
    	{
    		item->setCheckState(Qt::Checked);
    	}
    	//highlight active item
    	if (i == cc->activeLayerIndex())
    	{
		    layer_manager_->blockSignals(true);
				layer_manager_->setCurrentItem(item);
				layer_manager_->blockSignals(false);
    	}
    }
  }

	void TOPPViewBase::layerSelectionChange(int i)
	{
		if (i!=-1) activeCanvas_()->activateLayer(i);
	}

	void TOPPViewBase::layerContextMenu(const QPoint & pos)
	{
		QListWidgetItem* item = layer_manager_->itemAt(pos);
		if (item && item!=layer_manager_->item(0))
		{
			QMenu* context_menu = new QMenu(layer_manager_);
			context_menu->addAction("Delete");
			if (context_menu->exec(layer_manager_->mapToGlobal(pos)))
			{
				activeCanvas_()->removeLayer(layer_manager_->row(item));
				updateLayerbar();
			}
			delete (context_menu);
		}
	}

	void TOPPViewBase::layerVisibilityChange(QListWidgetItem* item)
	{
		int layer = layer_manager_->row(item);
		bool visible = activeCanvas_()->getLayer(layer).visible;
		
		if (item->checkState()==Qt::Unchecked && visible)
		{
			activeCanvas_()->changeVisibility(layer,false);
		}
		else if (item->checkState()==Qt::Checked && !visible)
		{
			activeCanvas_()->changeVisibility(layer,true);
		}
			
			
	}

  void TOPPViewBase::updateTabBar(QWidget* w)
  {
  	if (w)
  	{
  		tab_bar_->setCurrentIndex(dynamic_cast<SpectrumWidget*>(w)->window_id);
  	}
  }

  void TOPPViewBase::tileVertical()
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
      int actHeight = std::max(heightForEach, preferredHeight);

      window->parentWidget()->setGeometry( 0, y, ws_->width(), actHeight );
      y += actHeight;
    }
  }

  void TOPPViewBase::tileHorizontal()
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
      if ( window->windowState() & Qt::WindowMaximized )
      {
        // prevent flicker
        window->hide();
        window->showNormal();
      }
      int preferredWidth = window->minimumWidth()+window->parentWidget()->baseSize().width();

      int actWidth = std::max(widthForEach, preferredWidth);

      window->parentWidget()->setGeometry( y, 0, actWidth , ws_->height() );
      y += actWidth;
    }
  }

  void TOPPViewBase::showAsWindow_(SpectrumWidget* sw, const String& caption)
  {
  	ws_->addWindow(sw);
    connect(sw->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(updateToolbar()));
    connect(sw,SIGNAL(sendStatusMessage(std::string,OpenMS::UInt)),this,SLOT(showStatusMessage(std::string,OpenMS::UInt)));
    connect(sw,SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
    connect(sw,SIGNAL(modesChanged(QWidget*)),this,SLOT(updateToolbar()));
  
  	Spectrum2DWidget* sw2 = dynamic_cast<Spectrum2DWidget*>(sw);
  	if (sw2 != 0)
  	{
  		connect(sw2->getHorizontalProjection(),SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
  		connect(sw2->getVerticalProjection(),SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
  		connect(sw2,SIGNAL(showCurrentPeaksAs3D()),this,SLOT(showCurrentPeaksAs3D()));
  		connect(sw2,SIGNAL(showSpectrumAs1D(int)),this,SLOT(showSpectrumAs1D(int)));
  	}
  	
	  sw->setWindowTitle(caption.c_str());
    addTab_(sw,caption);
  }

  void TOPPViewBase::gotoDialog()
  {
    SpectrumWidget* w = activeWindow_();
    if (w)
    {
      w->showGoToDialog();
    }
  }

  SpectrumWidget*  TOPPViewBase::activeWindow_() const
  {
    return dynamic_cast<SpectrumWidget*>(ws_->activeWindow());
  }

  SpectrumCanvas*  TOPPViewBase::activeCanvas_() const
  {
    SpectrumWidget* sw = dynamic_cast<SpectrumWidget*>(ws_->activeWindow());
    if (sw == 0)
    {
    	return 0;
    }
    return sw->canvas();
  }

  Spectrum1DWidget* TOPPViewBase::active1DWindow_() const
  {
    Spectrum1DWidget* s1;
    if ((s1 = dynamic_cast<Spectrum1DWidget*>(ws_->activeWindow())))
    {
      return s1;
    }
    return 0;
  }

  Spectrum2DWidget* TOPPViewBase::active2DWindow_() const
  {
    Spectrum2DWidget* s2;
    if ((s2 = dynamic_cast<Spectrum2DWidget*>(ws_->activeWindow())))
    {
      return s2;
    }
    return 0;
  }

  Spectrum3DWidget* TOPPViewBase::active3DWindow_() const
  {
    Spectrum3DWidget* s3;
    if ((s3 = dynamic_cast<Spectrum3DWidget*>(ws_->activeWindow())))
    {
      return s3;
    }
    return 0;
  }

  void TOPPViewBase::loadPreferences(String filename)
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

    //load preferences, if file exists
    if (File::exists(filename))
    {
    	Param tmp;
    	tmp.load(filename);
    	//apply preferences if they are of the current TOPPView version
    	if(tmp.getValue("preferences:version").toString()==VersionInfo::getVersion())
    	{
      	setParameters(tmp);
    	}
    	else
    	{
    		cout << "The preferences files '" << filename  
    		     << "' is replaced as it is not of the current TOPPView version." << endl;
    	}
    }
    else
    {
      if (filename != default_ini_file)
      {
        cerr << "Unable to load INI File: '" << filename << "'" << endl;
      }
    }
    param_.setValue("PreferencesFile" , filename);

    //set the recent files
    Param p = param_.copy("preferences:RecentFiles");
    if (p.size()!=0)
    {
      for (Param::ParamIterator it=p.begin() ; it!=p.end() ; ++it)
      {
      	QString filename = it->value.toQString();
      	if (File::exists(filename.toAscii().data())) recent_files_.append(filename);
      }
    }
		
    updateRecentMenu_();
  }

  void TOPPViewBase::savePreferences()
  {
    // replace recent files
    param_.remove("preferences:RecentFiles");

    for (int i = 0; i < recent_files_.size(); ++i)
    {
      param_.setValue("preferences:RecentFiles:"+String(i),recent_files_[i]);
    }
		
		//set version
		param_.setValue("preferences:version",VersionInfo::getVersion());
		
    //save only the subsection that begins with "preferences:"
    try
    {
      param_.copy("preferences:").store(string(param_.getValue("PreferencesFile")));
    }
    catch(Exception::UnableToCreateFile& e)
    {
      cerr << "Unable to create INI File: '" << string(param_.getValue("PreferencesFile")) << "'" << endl;
    }
  }

  void TOPPViewBase::openRecentFile()
  {
		QAction* action = qobject_cast<QAction *>(sender());
    if (action)
		{
	  	setCursor(Qt::WaitCursor);
 	  	OpenDialog::Mower mow = OpenDialog::NO_MOWER;
			if ( (String)param_.getValue("preferences:intensity_cutoff")=="noise_estimator")
			{
				mow = OpenDialog::NOISE_ESTIMATOR;
			}
   		addSpectrum(action->text(),!recent_as_new_layer_->isChecked(),(String)param_.getValue("preferences:default_map_view")=="2d",true,mow);
			setCursor(Qt::ArrowCursor); 	
		}
 }

  void TOPPViewBase::openSpectrumDialog()
  {
    OpenDialog dialog(param_,this);
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

	void TOPPViewBase::editParamDialog()
	{
		// CREATE DIALOG
		QDialog dialog(this);
		QGridLayout* layout = new QGridLayout(&dialog);
		//Editor
		ParamEditor* edit = new ParamEditor(&dialog);
		edit->createShortcuts();
		layout->addWidget(edit,0,0,1,3);
		//Stretch
		layout->setColumnStretch(0,2);
		//Store button
		QPushButton* button = new QPushButton("Cancel",&dialog);
		connect(button,SIGNAL(pressed()),&dialog,SLOT(reject()));
		layout->addWidget(button,1,1);
		//Cancel button
		button = new QPushButton("OK",&dialog);
		connect(button,SIGNAL(pressed()),&dialog,SLOT(accept()));
		layout->addWidget(button,1,2);
		
		//LOAD DATA	
		QString name = QFileDialog::getOpenFileName(this,"Select a INI file",param_.getValue("preferences:default_path").toString().c_str(),tr("ini files (*.ini);; all files (*.*)"));
		if (name=="") return;

		Param p;
		p.load(name.toAscii().data());
		edit->loadEditable(p);

		//EXECUTE DIALOG + STORE DATA
		if (dialog.exec())
		{
			if (!File::writable(name.toAscii().data()))
			{
				QMessageBox::critical(this,"Error writing file!",(String("Cannot write to '")+name.toAscii().data()+"'!").c_str());
				return;
			}
			edit->store();
			p.store(name.toAscii().data());
		}
	}

	void TOPPViewBase::showTOPPDialog()
	{
		//check if there is a active window
		if (ws_->activeWindow())
		{
			String tmp_dir = param_.getValue("preferences:tmp_file_path").toString();
			String default_dir = param_.getValue("preferences:default_path").toString();
			
			ToolsDialog dialog(this,tmp_dir,default_dir,getCurrentLayer());
		
			if(dialog.exec()==QDialog::Accepted)
			{
				const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
				//warn if hidden layer => wrong layer selected...
				if (!layer.visible)
				{
					QMessageBox::warning(this,"Warning","The current layer is not visible!");
				}
				if (!File::writable(tmp_dir+"/in"))
				{
					QMessageBox::critical(this,"Error creating temporary file!",(String("Cannot write to '")+tmp_dir+"/in'!").c_str());
					return;
				}
				if (!File::writable(tmp_dir+"/out"))
				{
					QMessageBox::critical(this,"Error creating temporary file!",(String("Cannot write to '")+tmp_dir+"/out'!").c_str());
					return;
				}
				if (layer.type==LayerData::DT_PEAK)
				{
					MzDataFile().store(tmp_dir+"/in",layer.peaks);
				}
				else if (layer.type==LayerData::DT_FEATURE)
				{
					FeatureXMLFile().store(tmp_dir+"/in",layer.features);
				}
				else if (layer.type==LayerData::DT_FEATURE_PAIR)
				{
					LayerData::FeatureMapType feature_map=layer.features;
					FeatureMap<>::ConstIterator end=--feature_map.end();
					vector< ElementPair< Feature > > feature_pairs;
					FeatureMap<>::ConstIterator next;
					for(FeatureMap<>::ConstIterator i=feature_map.begin();i<end;i+=2)
					{
						next=++i;
						--i;
						feature_pairs.push_back(ElementPair<>(*i,*next));
					}
					FeaturePairsXMLFile().store(tmp_dir+"/in",feature_pairs);
				}
				else
				{
					return;
				}
				
				//compose argument list
				QStringList args;
				args <<"-ini" 
				        << QString("%1/in.ini").arg(tmp_dir.c_str())
				        << QString("-%1").arg(dialog.getInput().c_str())
				        <<  QString("%1/in").arg(tmp_dir.c_str());
				if (!dialog.noOutputAction())
				{
					args << QString("-%1").arg(dialog.getOutput().c_str())
				           <<  QString("%1/out").arg(tmp_dir.c_str());
				}
				//delete log and show it
				qobject_cast<QWidget *>(log_->parent())->show();
				log_->clear();
				//start process
				QProcess process;
				process.setProcessChannelMode(QProcess::MergedChannels);
				connect(&process,SIGNAL(readyReadStandardOutput()),this,SLOT(updateProcessLog()));
				process.start(dialog.getTool().c_str(),args);
				if (!process.waitForFinished(80000000))
				{
					QMessageBox::critical(this,"Execution of TOPP tool not successful!","The tool returned a exit code other than 0.<br>See log window for details.");
				}
				else if (process.exitStatus()==QProcess::CrashExit)
				{
					QMessageBox::critical(this,"Execution of TOPP tool not successful!",(String("The tool crashed during execution.<br>If you want to debug this crash, check the input files in '") + tmp_dir + "' or enable 'debug' mode in the TOPP ini file.").toQString());
				}
				else if (!File::readable(tmp_dir+"/out"))
				{
					QMessageBox::critical(this,"Error creating temporary file!",(String("Cannot read '")+tmp_dir+"/in'!").c_str());
					return;
				}
				else if(dialog.openAsWindow())
				{
					addSpectrum(tmp_dir+"/out",true,true,true);
				}
				else if(dialog.openAsLayer())
				{
					addSpectrum(tmp_dir+"/out",false,true,true);
				}
			}
		}
	}
	
	const LayerData* TOPPViewBase::getCurrentLayer() const
	{
		SpectrumCanvas* canvas = activeCanvas_();
		if (canvas==0)
		{
			return 0;
		}
		return &(canvas->getCurrentLayer());
	}

  void TOPPViewBase::showProjections()
  {
    Spectrum2DWidget* w = active2DWindow_();
    if (w)
    {
    	w->canvas()->showProjections();
    }
  }

	void TOPPViewBase::annotateWithID()
	{
		//check if there is a active window
		if (ws_->activeWindow())
		{
			const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
			//warn if hidden layer => wrong layer selected...
			if (!layer.visible)
			{
				QMessageBox::warning(this,"Warning","The current layer is not visible!");
			}
					
			//load id data
			QString name = QFileDialog::getOpenFileName(this,"Select ProteinIdentification data",param_.getValue("preferences:default_path").toString().c_str(),tr("identfication files (*.idXML);; all files (*.*)"));
			if(name!="")
			{
				vector<PeptideIdentification> identifications; 
				vector<ProteinIdentification> protein_identifications; 
				IdXMLFile().load(name, protein_identifications, identifications);
				if (layer.type==LayerData::DT_PEAK)
				{
					IDSpectrumMapper().annotate(const_cast<LayerData&>(layer).peaks, identifications, 0.1);
				}
				else if (layer.type==LayerData::DT_FEATURE || layer.type==LayerData::DT_FEATURE_PAIR)
				{
					IDFeatureMapper().annotate(const_cast<LayerData&>(layer).features,identifications,protein_identifications);
				}
			}
		}
	}

	void TOPPViewBase::showCurrentPeaksAs3D()
	{
  	//check if there is a active window
    if (ws_->activeWindow())
    {
      const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
    	const SpectrumCanvas::AreaType& area = activeWindow_()->canvas()->getVisibleArea();
    	const LayerData::ExperimentType& peaks = activeWindow_()->canvas()->getCurrentPeakData();
    	
    	if (layer.type==LayerData::DT_PEAK)
    	{
    		//open new 3D widget
    		Spectrum3DWidget* w = new Spectrum3DWidget(param_.copy("preferences:3d:",true), ws_);
  			SpectrumCanvas::ExperimentType& out = w->canvas()->addEmptyPeakLayer();
  			
    		for (LayerData::ExperimentType::ConstIterator it=peaks.RTBegin(area.min()[1]); it!=peaks.RTEnd(area.max()[1]); ++it)
    		{
    			if (it->getMSLevel()!=1) continue;
    			SpectrumCanvas::ExperimentType::SpectrumType spectrum;
    				
  				spectrum.setRT(it->getRT());
  				spectrum.setMSLevel(it->getMSLevel());
  				for (LayerData::ExperimentType::SpectrumType::ConstIterator it2 = it->MZBegin(area.min()[0]); it2!= it->MZEnd(area.max()[0]); ++it2)
  				{
  					if ( it2->getIntensity() >= layer.min_int && it2->getIntensity() <= layer.max_int)
  					{
  						spectrum.push_back(*it2);
  					}
  				}
  				out.push_back(spectrum);
    		}
    		
    		if (out.size()==0) //no extracted data
    		{
    		  QMessageBox::critical(this,"Error","The displayed region of the current layer is empty!");
    		  delete(w);
    		}
    		else //finish adding
    		{
    			String caption = layer.name + " (3D)";
    			w->canvas()->finishAdding(0.0);
					w->canvas()->setCurrentLayerName(caption);
		      showAsWindow_(w,caption);
		      w->showMaximized();
			    updateLayerbar();
    		}
    		
			}
    }
	}

	void TOPPViewBase::showSpectrumAs1D(int index)
	{
  	//check if there is a active window
    if (ws_->activeWindow())
    {
      const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
    	const LayerData::ExperimentType& peaks = activeWindow_()->canvas()->getCurrentPeakData();
    		
    	if (layer.type==LayerData::DT_PEAK)
    	{
    		//open new 1D widget
    		Spectrum1DWidget* w = new Spectrum1DWidget(param_.copy("preferences:1d:",true), ws_);

  			w->canvas()->addEmptyPeakLayer().push_back(peaks[index]);
  			String caption = layer.name + " (RT: " + String(peaks[index].getRT()) + ")";
  			w->canvas()->finishAdding(0.0);
				w->canvas()->setCurrentLayerName(caption);
	      showAsWindow_(w,caption);
	      w->showMaximized();
		    updateLayerbar();
			}
    }
	}

	void TOPPViewBase::showAboutDialog()
	{
		//dialog and grid layout 
		QDialog* dlg = new QDialog(this);
		QGridLayout* grid = new QGridLayout(dlg);
		dlg->setWindowTitle("About TOPPView");
		
		//image
		QLabel* label = new QLabel(dlg);
		QPixmap image(Oesterberg);
		label->setPixmap(image);
		grid->addWidget(label,0,0);
		
		//text
		QString text = QString("<BR>"
									 				 "<FONT size=+3>TOPPView</font><BR>"
									 				 "<BR>"
													 "Version: %1<BR>"
													 "<BR>"
													 "OpenMS and TOPP is free software available under the<BR>"
													 "Lesser GNU Public License (LGPL)").arg(VersionInfo::getVersion().c_str());
		label = new QLabel(text,dlg);
		grid->addWidget(label,0,1,Qt::AlignTop | Qt::AlignLeft);
		
		//execute
		dlg->exec();
	}

	void TOPPViewBase::renameLayer()
	{
		//check if there is a active window
		if (ws_->activeWindow())
		{
			const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
			//warn if hidden layer => wrong layer selected...
			if (!layer.visible)
			{
				QMessageBox::warning(this,"Warning","The current layer is not visible!");
			}
			QString name = QInputDialog::getText(this,"New layer name","Name:");
			if (name!="")
			{
				const_cast<LayerData&>(layer).name = name;
				updateLayerbar();
			}
		}
	}

	void TOPPViewBase::updateProcessLog()
	{
		QProcess* process = qobject_cast<QProcess *>(sender());
		log_->textCursor().insertText(process->readAllStandardOutput());
	}

} //namespace OpenMS

