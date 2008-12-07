// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/VISUAL/DIALOGS/DataFilterDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/DBOpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>
#include <OpenMS/VISUAL/DIALOGS/SpectrumAlignmentDialog.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum3DWidget.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/BASELINE/TopHatFilter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewPrefDialog.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>

//Qt
#include <QtGui/QToolBar>
#include <QtGui/QDockWidget>
#include <QtGui/QListWidget>
#include <QtGui/QListWidgetItem>
#include <QtGui/QTreeWidget>
#include <QtGui/QTreeWidgetItem>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QToolButton>
#include <QtGui/QMessageBox>
#include <QtGui/QListWidgetItem>
#include <QtGui/QToolTip>
#include <QtGui/QFileDialog>
#include <QtGui/QPainter>
#include <QtCore/QDir>
#include <QtCore/QDate>
#include <QtGui/QWhatsThis>
#include <QtGui/QInputDialog>
#include <QtGui/QTextEdit>
#include <QtGui/QCheckBox>
#include <QtGui/QCloseEvent>
#include <QtGui/QDesktopServices>
#include <QtCore/QUrl>
#include <QtGui/QSplashScreen>

//intensity modes
#include "../VISUAL/ICONS/lin.xpm"
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
#include "../VISUAL/ICONS/convexhull.xpm"
#include "../VISUAL/ICONS/convexhulls.xpm"
#include "../VISUAL/ICONS/numbers.xpm"
#include "../VISUAL/ICONS/elements.xpm"

//misc
#include "../VISUAL/ICONS/TOPPView.xpm"
#include "../VISUAL/ICONS/Oesterberg.xpm"

#include <algorithm>
#include <utility>

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
    tab_bar_->setWhatsThis("Tab bar<BR><BR>Close tabs through the context menu or by double-clicking them.<BR>The tab bar accepts drag-and-drop from the layer bar.");
    tab_bar_->addTab("dummy",4710);
    tab_bar_->setMinimumSize(tab_bar_->sizeHint());
    tab_bar_->removeId(4710);
    //connect slots and sigals for selecting spectra
    connect(tab_bar_,SIGNAL(currentIdChanged(int)),this,SLOT(focusByTab(int)));
    connect(tab_bar_,SIGNAL(aboutToCloseId(int)),this,SLOT(closeByTab(int)));
		//connect signals ans slots for drag-and-drop
		connect(tab_bar_,SIGNAL(dropOnWidget(const QMimeData*)),this,SLOT(copyLayer(const QMimeData*)));		
		connect(tab_bar_,SIGNAL(dropOnTab(const QMimeData*,int)),this,SLOT(copyLayer(const QMimeData*, int)));		

    box_layout->addWidget(tab_bar_);
    ws_=new QWorkspace(dummy);
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateToolBar()));
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateTabBar(QWidget*)));
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateLayerBar()));
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateSpectrumBar()));
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateFilterBar()));
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateMenu()));
 
    box_layout->addWidget(ws_);

		//################## MENUS #################
    // File menu
    QMenu* file = new QMenu("&File",this);
    menuBar()->addMenu(file);
    file->addAction("&Open file",this,SLOT(openFileDialog()), Qt::CTRL+Qt::Key_O);
    file->addAction("&Open from database",this,SLOT(openDatabaseDialog()), Qt::CTRL+Qt::Key_D);
    file->addAction("&Close",this,SLOT(closeFile()), Qt::CTRL+Qt::Key_W);
		file->addSeparator();
		
		//Meta data
		file->addAction("&Show meta data (file)",this,SLOT(metadataFileDialog()));
    file->addAction("&Show meta data (database)",this,SLOT(metadataDatabaseDialog()));
    file->addSeparator();
		
		//Recent files    
    QMenu* recent_menu = new QMenu("&Recent files", this);
  	recent_actions_.resize(20);
		for (UInt i = 0; i<20; ++i)
		{
			recent_actions_[i] = recent_menu->addAction("",this,SLOT(openRecentFile()));
			recent_actions_[i]->setVisible(false);
		}
  	file->addMenu(recent_menu);

    file->addSeparator();
    file->addAction("&Preferences",this, SLOT(preferencesDialog()));
    file->addAction("&Quit",qApp,SLOT(quit()));
    
    //Tools menu
    QMenu* tools = new QMenu("&Tools",this);
    menuBar()->addMenu(tools);
    tools->addAction("&Go to",this,SLOT(gotoDialog()), Qt::CTRL+Qt::Key_G);
    tools->addAction("&Edit meta data",this,SLOT(editMetadata()), Qt::CTRL+Qt::Key_M);
    tools->addAction("&Statistics",this,SLOT(layerStatistics()));
		tools->addSeparator();
    tools->addAction("Apply TOPP tool (whole layer)", this, SLOT(showTOPPDialog()), Qt::CTRL+Qt::Key_T)->setData(false);
    tools->addAction("Apply TOPP tool (visible layer data)", this, SLOT(showTOPPDialog()), Qt::CTRL+Qt::SHIFT+Qt::Key_T)->setData(true);
    tools->addAction("Rerun TOPP tool", this, SLOT(rerunTOPPTool()),Qt::Key_F4);
    tools->addSeparator();
    tools->addAction("&Annotate with identifiction", this, SLOT(annotateWithID()), Qt::CTRL+Qt::Key_I);
    tools->addAction("Align spectra", this, SLOT(showSpectrumAlignmentDialog()));
    tools->addAction("Generate theoretical spectrum", this, SLOT(showSpectrumGenerationDialog()));

    //Layer menu
    QMenu* layer = new QMenu("&Layer",this);
    menuBar()->addMenu(layer);
    layer->addAction("Save all data", this, SLOT(saveLayerAll()), Qt::CTRL+Qt::Key_S);
    layer->addAction("Save visible data", this, SLOT(saveLayerVisible()), Qt::CTRL+Qt::SHIFT+Qt::Key_S);
		layer->addSeparator();
    layer->addAction("Show/hide grid lines", this, SLOT(toggleGridLines()), Qt::CTRL+Qt::Key_R);
    layer->addAction("Show/hide axis legends", this, SLOT(toggleAxisLegends()), Qt::CTRL+Qt::Key_L);
		layer->addSeparator();
    layer->addAction("Preferences", this, SLOT(showPreferences()));
    
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
		QAction* action = help->addAction("OpenMS website",this,SLOT(showURL()));
		action->setData("http://www.OpenMS.de");
		action = help->addAction("TOPPView tutorial (online)",this,SLOT(showURL()), Qt::Key_F1);
		action->setData("http://www-bs2.informatik.uni-tuebingen.de/services/OpenMS-release/html/TOPPViewTutorial.html");		
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
    
    //intensity modes
    intensity_group_ = new QButtonGroup(tool_bar_);
    intensity_group_->setExclusive(true);
    
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(lin));
    b->setToolTip("Intensity: Normal");
    b->setShortcut(Qt::Key_N);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Normal<BR><BR>Intensity is displayed unmodified.<BR>(Hotkey: N)");
    intensity_group_->addButton(b,SpectrumCanvas::IM_NONE);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(percentage));
    b->setToolTip("Intensity: Percentage");
    b->setShortcut(Qt::Key_P);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Percentage<BR><BR>Intensity is displayed as a percentage of the layer"
    								" maximum intensity. If only one layer is displayed this mode behaves like the"
    								" normal mode. If more than one layer is displayed intensities are aligned."
    								"<BR>(Hotkey: P)");
    intensity_group_->addButton(b,SpectrumCanvas::IM_PERCENTAGE);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(snap));
    b->setToolTip("Intensity: Snap to maximum displayed intensity");
    b->setShortcut(Qt::Key_S);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Snap to maximum displayed intensity<BR><BR> In this mode the"
    								" color gradient is adapted to the maximum currently displayed intensity."
    								"<BR>(Hotkey: S)");
    intensity_group_->addButton(b,SpectrumCanvas::IM_SNAP);
		tool_bar_->addWidget(b);
    connect(intensity_group_,SIGNAL(buttonClicked(int)),this,SLOT(setIntensityMode(int)));
    tool_bar_->addSeparator();

    //common buttons
    QAction* reset_zoom_button = tool_bar_->addAction(QPixmap(reset_zoom), "Reset Zoom", this, SLOT(resetZoom()));
    reset_zoom_button->setWhatsThis("Reset zoom: Zooms out as far as possible and resets the zoom history.<BR>(Hotkey: Backspace)");

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

    dm_precursors_2d_ = tool_bar_2d_->addAction(QPixmap(precursors),"Show fragment scan precursors");
    dm_precursors_2d_->setCheckable(true);
    dm_precursors_2d_->setWhatsThis("2D peak draw mode: Precursors<BR><BR>fragment scan precursor peaks are marked");
    connect(dm_precursors_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    projections_2d_ = tool_bar_2d_->addAction(QPixmap(projections), "Show Projections" ,this, SLOT(toggleProjections()));
    projections_2d_->setWhatsThis("Projections: Shows projections of peak data along RT and MZ axis.");

    dm_hull_2d_ = tool_bar_2d_->addAction(QPixmap(convexhull),"Show feature convex hull");
    dm_hull_2d_->setCheckable(true);
    dm_hull_2d_->setWhatsThis("2D feature draw mode: Convex hull<BR><BR>The convex hull of the feature is displayed");
    connect(dm_hull_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    dm_hulls_2d_ = tool_bar_2d_->addAction(QPixmap(convexhulls),"Show feature convex hulls");
    dm_hulls_2d_->setCheckable(true);
    dm_hulls_2d_->setWhatsThis("2D feature draw mode: Convex hulls<BR><BR>The convex hulls of the feature are displayed: One for each mass trace.");
    connect(dm_hulls_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    dm_numbers_2d_ = tool_bar_2d_->addAction(QPixmap(numbers),"Show feature identifiers");
    dm_numbers_2d_->setCheckable(true);
    dm_numbers_2d_->setWhatsThis("2D feature draw mode: Numbers/labels<BR><BR>The feature number is displayed next to the feature. If the meta data value 'label' is set, it is displayed in brackets after the number.");
    connect(dm_numbers_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    dm_elements_2d_ = tool_bar_2d_->addAction(QPixmap(elements),"Show consensus feature element positions");
    dm_elements_2d_->setCheckable(true);
    dm_elements_2d_->setWhatsThis("2D consensus feature draw mode: Elements<BR><BR>The individual elements that make up the  consensus feature are drawn.");
    connect(dm_elements_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));


		//################## Dock widgets #################
    //layer window
    QDockWidget* layer_bar = new QDockWidget("Layers", this);
    addDockWidget(Qt::RightDockWidgetArea, layer_bar);
    layer_manager_ = new QListWidget(layer_bar);
    layer_manager_->setWhatsThis("Layer bar<BR><BR>Here the availabe layers are shown. Left-click on a layer to select it.<BR>Layers can be shown and hidden using the checkboxes in front of the name.<BR> Renaming and removing a layer is possible through the context menu.<BR>Dragging a layer to the tab bar copies the layer.<BR>Double-clicking a layer open its preferences.<BR>You can use the 'PageUp' and 'PageDown' buttons to change the selected layer.");

    layer_bar->setWidget(layer_manager_);
    layer_manager_->setContextMenuPolicy(Qt::CustomContextMenu);
		layer_manager_->setDragEnabled(true);
    connect(layer_manager_,SIGNAL(currentRowChanged(int)),this,SLOT(layerSelectionChange(int)));
		connect(layer_manager_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(layerContextMenu(const QPoint&)));
		connect(layer_manager_,SIGNAL(itemChanged(QListWidgetItem*)),this,SLOT(layerVisibilityChange(QListWidgetItem*)));
    connect(layer_manager_,SIGNAL(itemDoubleClicked(QListWidgetItem*)),this,SLOT(layerEdit(QListWidgetItem*)));
    
    windows->addAction("&Show layer window",layer_bar,SLOT(show()));
		
    //spectrum selection
    QDockWidget* spectrum_bar = new QDockWidget("Spectra", this);
    addDockWidget(Qt::RightDockWidgetArea, spectrum_bar);
    spectrum_selection_ = new QTreeWidget(spectrum_bar);
    spectrum_selection_->setWhatsThis("Spectrum selection bar<BR><BR>Here all spectra of the current experiment are shown. Left-click on a spectrum to open it.");
    spectrum_selection_->setColumnCount(3);
  	QStringList header_labels;
  	header_labels.append(QString("MS level"));
  	header_labels.append(QString("RT"));
  	header_labels.append(QString("m/z"));
  	spectrum_selection_->setHeaderLabels(header_labels);
    spectrum_bar->setWidget(spectrum_selection_);
    connect(spectrum_selection_,SIGNAL(itemClicked(QTreeWidgetItem*, int)),this,SLOT(spectrumSelectionChange(QTreeWidgetItem*, int)));
    
    windows->addAction("&Show spectrum selection window",spectrum_bar,SLOT(show()));
        
    //data filters
    QDockWidget* filter_bar = new QDockWidget("Data filters", this);
    addDockWidget(Qt::RightDockWidgetArea, filter_bar);
    QWidget* tmp_widget = new QWidget(); //dummy widget as QDockWidget takes only one widget
    filter_bar->setWidget(tmp_widget);
    
    QVBoxLayout* vbl = new QVBoxLayout(tmp_widget);
      
    filters_ = new QListWidget(tmp_widget);
    filters_->setSelectionMode(QAbstractItemView::NoSelection);
    filters_->setWhatsThis("Data filter bar<BR><BR>Here filtering options for the current layer can be set.<BR>Through the context menu you can add, remove and edit filters.<BR>For convenience, editing filters is also possible by double-clicking them.");
    filters_->setContextMenuPolicy(Qt::CustomContextMenu);
		connect(filters_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(filterContextMenu(const QPoint&)));
		connect(filters_,SIGNAL(itemDoubleClicked(QListWidgetItem*)),this,SLOT(filterEdit(QListWidgetItem*)));
		vbl->addWidget(filters_);

    filters_check_box_ = new QCheckBox("Enable/disable all filters", tmp_widget);
    connect(filters_check_box_,SIGNAL(toggled(bool)),this,SLOT(layerFilterVisibilityChange(bool)));
    vbl->addWidget(filters_check_box_); 
    
    windows->addAction("&Show filter window",filter_bar,SLOT(show()));


		//log window
		QDockWidget* log_bar = new QDockWidget("Log", this);
		addDockWidget(Qt::BottomDockWidgetArea, log_bar);
		log_ = new QTextEdit(log_bar);
		log_->setReadOnly(true);
		log_->setContextMenuPolicy(Qt::CustomContextMenu);
		connect(log_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(logContextMenu(const QPoint&)));
		log_bar->setWidget(log_);
		log_bar->hide();
    windows->addAction("&Show log window",log_bar,SLOT(show()));
		
		//################## DEFAULTS #################
    //general
    defaults_.setValue("preferences:default_map_view", "2d", "Default visualization mode for maps.");
		defaults_.setValidStrings("preferences:default_map_view",StringList::create("2d,3d"));
    defaults_.setValue("preferences:default_path", ".", "Default path for loading and storing files.");
    defaults_.setValue("preferences:tmp_file_path", QDir::tempPath(), "Path where temporary files can be created.");
    defaults_.setValue("preferences:number_of_recent_files", 15, "Number of recent files in the main menu.");
    defaults_.setMinInt("preferences:number_of_recent_files",5);
    defaults_.setMaxInt("preferences:number_of_recent_files",20);
    defaults_.setValue("preferences:legend", "show", "Legend visibility");
		defaults_.setValidStrings("preferences:legend",StringList::create("show,hide"));
    defaults_.setValue("preferences:intensity_cutoff", "off","Low intensity cutoff for maps.");
		defaults_.setValidStrings("preferences:intensity_cutoff",StringList::create("on,off"));
    defaults_.setValue("preferences:on_file_change","ask","What action to take, when a data file changes. Do nothing, update automatically or ask the user.");
		defaults_.setValidStrings("preferences:on_file_change",StringList::create("none,ask,update automatically"));
    //db
    defaults_.setValue("preferences:db:host", "localhost", "Database server host name.");
    defaults_.setValue("preferences:db:login", "NoName", "Database login.");
    defaults_.setValue("preferences:db:name", "OpenMS", "Database name.");
    defaults_.setValue("preferences:db:port", 3306, "Database server port.");
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
  	
  	//update the menu
  	updateMenu();
		
		//######################### Additional context menus ######################################
		add_2d_context_ = new QMenu("More",this);
		add_2d_context_->addAction("Show layer in 3D",this,SLOT(showCurrentPeaksAs3D()));
  
		topp_.process = 0;
	}

  TOPPViewBase::~TOPPViewBase()
  {
  	savePreferences();
  	abortTOPPTool();
  }

  void TOPPViewBase::closeEvent(QCloseEvent* event)
  {
  	ws_->closeAllWindows();
  	event->accept();
  }

	void TOPPViewBase::showURL()
	{
		QAction* action = qobject_cast<QAction*>(sender());
		if (!QDesktopServices::openUrl(QUrl(action->data().toString())))
		{
			QMessageBox::warning(this, tr("Error"),
												 tr("Unable to open\n") + 
														action->data().toString() + 
														tr("\n\nPossible reason: security settings or misconfigured Operating System"));
		}
	}

  void TOPPViewBase::addDataDB(UInt db_id, bool show_options, String caption, UInt window_id)
  {
  	//set wait cursor
  	setCursor(Qt::WaitCursor);
				
    //Open DB connection
    DBConnection con;
    connectToDB_(con);
    if (!con.isConnected())
    {
      setCursor(Qt::ArrowCursor);
      return;
    }

		//load the data
    DBAdapter db(con);
    ExperimentType exp;
    FeatureMapType dummy_map;
    ConsensusMapType dummy_map2;
		try
		{
			db.loadExperiment(db_id, exp);
		}
		catch (Exception::BaseException& e)
		{
			QMessageBox::critical(this,"Error",(String("Error while reading data: ")+e.what()).c_str());
      setCursor(Qt::ArrowCursor);
      return;
		}
		
		//determine if the data is 1D or 2D
    QSqlQuery result = con.executeQuery(String("SELECT count(id) from DATA_Spectrum where fid_MSExperiment='")+db_id+"' and MSLevel='1'");
		bool is_2D = (result.value(0).toInt()>1);
		
		//add data
    if (caption=="") caption = String("DB entry ")+db_id;
		addData_(dummy_map, dummy_map2, exp, false, is_2D, show_options, "", caption, window_id);
  	
  	//Reset cursor
  	setCursor(Qt::ArrowCursor);
  }

  float TOPPViewBase::estimateNoise_(const ExperimentType& exp)
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

      if (scan < exp.size() && exp[scan].getMSLevel()==1 && exp[scan].size()!=0)
      {
        vector<float> tmp;
        tmp.reserve(exp[scan].size());
        for(SpectrumType::ConstIterator it = exp[scan].begin()
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
    return noise / 10.0f;
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
		QComboBox* on_file_change = dlg.findChild<QComboBox*>("on_file_change");
		
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
		
		//set General values
		default_path->setText(param_.getValue("preferences:default_path").toQString());
		temp_path->setText(param_.getValue("preferences:tmp_file_path").toQString());
		recent_files->setValue((Int)param_.getValue("preferences:number_of_recent_files"));
		map_default->setCurrentIndex(map_default->findText(param_.getValue("preferences:default_map_view").toQString()));
		map_cutoff->setCurrentIndex(map_cutoff->findText(param_.getValue("preferences:intensity_cutoff").toQString()));		
		on_file_change->setCurrentIndex(on_file_change->findText(param_.getValue("preferences:on_file_change").toQString()));		


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
		
		//execute dialog
		if (dlg.exec())
		{
			param_.setValue("preferences:default_path", default_path->text());
			param_.setValue("preferences:tmp_file_path", temp_path->text());
			param_.setValue("preferences:number_of_recent_files", recent_files->value());
			param_.setValue("preferences:default_map_view", map_default->currentText());
			param_.setValue("preferences:intensity_cutoff", map_cutoff->currentText());
			param_.setValue("preferences:on_file_change", on_file_change->currentText());

			param_.setValue("preferences:db:host",db_host->text());
			param_.setValue("preferences:db:port",db_port->value());
			param_.setValue("preferences:db:name",db_name->text());
			param_.setValue("preferences:db:login",db_login->text());
			param_.remove("DBPassword");

			param_.setValue("preferences:1d:peak_color",color_1D->getColor().name());
			param_.setValue("preferences:1d:highlighted_peak_color",selected_1D->getColor().name());
			param_.setValue("preferences:1d:icon_color",icon_1D->getColor().name());

			param_.setValue("preferences:2d:dot:gradient",peak_2D->gradient().toString());
			param_.setValue("preferences:2d:mapping_of_mz_to",mapping_2D->currentText());

			param_.setValue("preferences:3d:dot:gradient",peak_3D->gradient().toString());
			param_.setValue("preferences:3d:dot:shade_mode", shade_3D->currentIndex());
			param_.setValue("preferences:3d:dot:line_width",line_width_3D->value());

			savePreferences();
		}
  }

  void TOPPViewBase::addDataFile(const String& filename,bool show_options, bool add_to_recent, String caption, UInt window_id)
  {
    setCursor(Qt::WaitCursor);

  	String abs_filename = File::absolutePath(filename);

  	//check if the file exists
    if (!File::exists(abs_filename))
    {
    	showLogMessage_(LS_ERROR,"Open file error",String("The file '")+abs_filename+"' does not exist!");
    	setCursor(Qt::ArrowCursor);
      return;
    }

		//determine file type
  	FileHandler fh;
		FileHandler::Type file_type = fh.getType(abs_filename);
		if (file_type==FileHandler::UNKNOWN)
		{
			showLogMessage_(LS_ERROR,"Open file error",String("Could not determine file type of '")+abs_filename+"'!");
    	setCursor(Qt::ArrowCursor);
      return;
		}
		//abort if file type unsupported
		if (file_type==FileHandler::PARAM || file_type==FileHandler::IDXML)
		{
			showLogMessage_(LS_ERROR,"Open file error",String("The type '")+fh.typeToName(file_type)+"' is not supported!");
   		setCursor(Qt::ArrowCursor);
      return;
		}
		
		//try to load data and determine if it is 1D or 2D data
		FeatureMapType feature_map;
		ExperimentType peak_map;
		ConsensusMapType consensus_map;
		
		bool is_2D = false;
		bool is_feature = false;
    try
    {
	    if (file_type==FileHandler::FEATUREXML)
	    {
        FeatureXMLFile().load(abs_filename,feature_map);
        is_2D = true;
        is_feature = true;
      }
      else if (file_type==FileHandler::CONSENSUSXML)
	    {
        ConsensusXMLFile().load(abs_filename,consensus_map);
        is_2D = true;
        is_feature = true;
      }
      else
      {
      	fh.loadExperiment(abs_filename,peak_map, file_type,ProgressLogger::GUI);
      	UInt ms1_scans = 0;
      	for (UInt i=0; i<peak_map.size();++i)
      	{
      		if (peak_map[i].getMSLevel()==1) ++ms1_scans;
      		if (ms1_scans>1)
      		{
      			is_2D = true;
      			break;
      		}
      	}
      }
    }
    catch(Exception::BaseException& e)
    {
    	showLogMessage_(LS_ERROR,"Error while loading file",e.what());
    	setCursor(Qt::ArrowCursor);
      return;
    }
    
    //try to add the data
		if (caption=="")
		{
			caption = File::basename(abs_filename);
    }
    else
    {
    	abs_filename = "";
    }
    addData_(feature_map, consensus_map, peak_map, is_feature, is_2D, show_options, abs_filename, caption, window_id);
  	
  	//add to recent file
  	if (add_to_recent) addRecentFile_(filename);
    
    //reset cursor
    setCursor(Qt::ArrowCursor);
  }
  
  void TOPPViewBase::addData_(FeatureMapType& feature_map, ConsensusMapType& consensus_map, ExperimentType& peak_map, bool is_feature, bool is_2D, bool show_options, const String& filename, const String& caption, UInt window_id)
  {
  	//initialize flags with defaults from the parameters
  	bool as_new_window = true;
  	bool maps_as_2d = ((String)param_.getValue("preferences:default_map_view")=="2d");
  	bool use_mower = ((String)param_.getValue("preferences:intensity_cutoff")=="on");
		
		bool is_consensus_feature = (feature_map==FeatureMapType());
		
		//set the window where (new layer) data could be opened in
		SpectrumWidget* open_window = window_(window_id);
		if (open_window==0)
		{
			open_window = activeWindow_();
		}
		else
		{
			as_new_window = false;
		}
		
		//create dialog no matter if it is shown or not. It is used to determine the flags
		TOPPViewOpenDialog dialog(caption, as_new_window, maps_as_2d, use_mower, this);
		//disable opening in new window when
		if (open_window==0 // there is no active window
			  || (is_2D && qobject_cast<Spectrum1DWidget*>(open_window)!=0) //2D data is to be opened, but the current window is a 1D window
			  || (is_feature && qobject_cast<Spectrum3DWidget*>(open_window)!=0)) //feature data is to be opened, but the current window is a 3D window
		{
			dialog.disableLocation(true);
		}
		//disable 2d/3d option for features and single scans
		if (is_feature || !is_2D) dialog.disableDimension(true);
		//disable cutoff for features and single scans
		if (is_feature || !is_2D) dialog.disableCutoff(false);
		//enable merge layers if a feature layer is opened and there are already features layers to merge it to
		if (is_feature && open_window!=0) //TODO merge
		{
			SpectrumCanvas* open_canvas = open_window->canvas();
			Map<UInt,String> layers;
			for (UInt i=0; i<open_canvas->getLayerCount(); ++i)
			{
				if (!is_consensus_feature && open_canvas->getLayer(i).type==LayerData::DT_FEATURE)
				{
					layers[i] = open_canvas->getLayer(i).name;
				}
				if (is_consensus_feature && open_canvas->getLayer(i).type==LayerData::DT_CONSENSUS)
				{
					layers[i] = open_canvas->getLayer(i).name;
				}
			}
			dialog.setMergeLayers(layers);
		}
		
		//show options if requested
		if (show_options && !dialog.exec())
		{
			return;
		}
		as_new_window = dialog.openAsNewWindow();
		maps_as_2d = dialog.viewMapAs2D();
		use_mower = dialog.isCutoffEnabled();
  	Int merge_layer = dialog.getMergeLayer();
  	
		//determine the window to open the data in
		if (as_new_window) //new window
    {
      if (!is_2D) //1d
      {
        open_window = new Spectrum1DWidget(getSpectrumParameters_(1), ws_);
      }
      else if (maps_as_2d || is_feature) //2d or features
      {
      	open_window = new Spectrum2DWidget(getSpectrumParameters_(2), ws_);
      	open_window->canvas()->setAdditionalContextMenu(add_2d_context_);
      } 
      else //3d
      {
      	open_window = new Spectrum3DWidget(getSpectrumParameters_(3), ws_);
      }
    }
		
    if (merge_layer==-1) //add data to the window
    {
	    if (is_feature) //features and consensus features
	    {
	    	if (!is_consensus_feature) //features
	    	{
	      	if (!open_window->canvas()->addLayer(feature_map,filename)) return;
	    	}
	    	else //consensus features
	    	{
	    		if (!open_window->canvas()->addLayer(consensus_map,filename)) return;
	    	}
	    }
	    else //peaks
	    {
			  if (!open_window->canvas()->addLayer(peak_map,filename)) return;
	      //calculate noise
	      if(use_mower && is_2D)
	      {
	        DoubleReal cutoff = estimateNoise_(open_window->canvas()->getCurrentLayer().peaks);
					//create filter
					DataFilters::DataFilter filter;
					filter.field = DataFilters::INTENSITY;
					filter.op = DataFilters::GREATER_EQUAL;
					filter.value = cutoff;
					///add filter
					DataFilters filters;
					filters.add(filter);
					open_window->canvas()->setFilters(filters);
	      }
			}
			
			//set caption
    	open_window->canvas()->setLayerName(open_window->canvas()->activeLayerIndex(), caption);
		}
		else //merge features data into feature layer
		{
			Spectrum2DCanvas* canvas = qobject_cast<Spectrum2DCanvas*>(open_window->canvas());
			if (is_consensus_feature)
			{
				canvas->mergeIntoLayer(merge_layer,consensus_map);
			}
			else
			{
				canvas->mergeIntoLayer(merge_layer,feature_map);
			}
		}

    if (as_new_window)
    {
    	showAsWindow_(open_window,caption);
		}
		updateLayerBar();
		updateSpectrumBar();
		updateFilterBar();
  	updateMenu();
	}
	
  void TOPPViewBase::addRecentFile_(const String& filename)
  {
  	//find out absolute path
    String tmp = File::absolutePath(filename);
	
		// remove the new file if already in the recent list and prepend it
		recent_files_.removeAll(tmp.toQString());
		recent_files_.prepend(tmp.toQString());
		
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

  SpectrumWidget* TOPPViewBase::window_(int id) const
  {
  	//cout << "Looking for tab with id: " << id << endl;
  	QList<QWidget*> windows = ws_->windowList();
		for(int i=0; i< windows.size(); ++i)
		{
			SpectrumWidget* window = qobject_cast<SpectrumWidget*>(windows.at(i));
			//cout << "  Tab " << i << ": " << window->window_id << endl;
			if (window->window_id == id)
			{
				return window;
			}
		}
		return 0;
  }

  void TOPPViewBase::closeByTab(int id)
  {
  	SpectrumWidget* window = window_(id);
  	if (window)
  	{
  		window->close();
  		updateMenu();
  	}
  }
 
  void TOPPViewBase::focusByTab(int id)
  {
  	SpectrumWidget* window = window_(id);
  	if (window)
  	{
  		window->setFocus();
  	}
  }
  
  void TOPPViewBase::closeFile()
  {
    ws_->activeWindow()->close();
    updateMenu();
  }

  void TOPPViewBase::editMetadata()
  {
  	SpectrumCanvas* canvas = activeCanvas_();
  	
    //warn if hidden layer => wrong layer selected...
  	if (!canvas->getCurrentLayer().visible)
  	{
  		showLogMessage_(LS_NOTICE,"The current layer is not visible","Have you selected the right layer for this action?");
  	}
  	
  	//show editable meta data dialog
  	canvas->showMetaData(true);
  }

  void TOPPViewBase::layerStatistics()
  {
  	activeWindow_()->showStatistics();
  	updateFilterBar();
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
    else if (isinf(mz) || isnan(mz))
		{
      int_label_->setText("m/z: n/a");
		}
    else
    {
      mz_label_->setText((String("m/z: ")+String::number(mz,3).fillLeft(' ',8)).toQString());
    }
    if (rt==-1)
    {
      rt_label_->setText("RT: ");
    }
    else if (isinf(rt) || isnan(rt))
		{
      int_label_->setText("RT: n/a");
		}
    else
    {
      rt_label_->setText((String("RT: ")+String::number(rt,1).fillLeft(' ',8)).toQString());
    }
    if (intensity==-1)
    {
      int_label_->setText("Int: ");
    }
    else if (isinf(intensity) || isnan(intensity))
		{
      int_label_->setText("Int: n/a");
		}
		else
    {
      int_label_->setText((String("Int: ")+String::number(intensity,1).fillLeft(' ',12)).toQString());
    }
    statusBar()->update();
  }

  void TOPPViewBase::resetZoom()
  {
    SpectrumWidget* window = activeWindow_();
    if (window!=0)
    {
      window->canvas()->resetZoom();
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
			else if (action == dm_hulls_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::F_HULLS,on);
			}
			else if (action == dm_hull_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::F_HULL,on);
			}
			else if (action == dm_numbers_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::F_NUMBERS,on);
			}
			//consensus features
			else if (action == dm_elements_2d_)
			{
				 win->canvas()->setLayerFlag(LayerData::C_ELEMENTS,on);
			}
		}
  }

  void TOPPViewBase::updateToolBar()
  {
    SpectrumWidget* w = activeWindow_();

    if (w)
    {
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
        Spectrum1DWidget* window = qobject_cast<Spectrum1DWidget*>(windows.at(i));
        if (window !=0 && window!=w)
        {
          link_box_->insertItem(++item_index,File::basename(window->windowTitle()).toQString(),window->window_id);
        	if (link_map_.find(w1->window_id)!=link_map_.end() && link_map_[w1->window_id] == window->window_id)
          {
            link_box_->setCurrentIndex(item_index);
          }
        }
      }

      //show/hide toolbars and buttons
      tool_bar_1d_->show();
      tool_bar_2d_->hide();
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
      	dm_hulls_2d_->setVisible(false);
      	dm_hull_2d_->setVisible(false);
      	dm_numbers_2d_->setVisible(false);
      	dm_elements_2d_->setVisible(false);
				dm_precursors_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::P_PRECURSORS));
			}
			//feature draw modes
			else if (w2->canvas()->getCurrentLayer().type == LayerData::DT_FEATURE)
			{
      	dm_precursors_2d_->setVisible(false);
      	projections_2d_->setVisible(false);
      	dm_hulls_2d_->setVisible(true);
      	dm_hull_2d_->setVisible(true);
      	dm_numbers_2d_->setVisible(true);
      	dm_elements_2d_->setVisible(false);
      	dm_hulls_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::F_HULLS));
      	dm_hull_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::F_HULL));
      	dm_numbers_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::F_NUMBERS));
			}
			//consensus feature draw modes
			else
			{
      	dm_precursors_2d_->setVisible(false);
      	projections_2d_->setVisible(false);
      	dm_hulls_2d_->setVisible(false);
      	dm_hull_2d_->setVisible(false);
      	dm_numbers_2d_->setVisible(false);
      	dm_elements_2d_->setVisible(true);
      	dm_hulls_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::C_ELEMENTS));
			}
      //show/hide toolbars and buttons
      tool_bar_1d_->hide();
      tool_bar_2d_->show();
    }

    //1D
    Spectrum3DWidget* w3 = active3DWindow_();
    if (w3)
    {
      //show/hide toolbars and buttons
      tool_bar_1d_->hide();
      tool_bar_2d_->hide();
    }
  }

  void TOPPViewBase::updateLayerBar()
  {
		layer_manager_->clear();
    SpectrumCanvas* cc = activeCanvas_();
    if (cc == 0)
    {
      return;
    }
		layer_manager_->blockSignals(true);
		QListWidgetItem* item = 0;
		QString name;
    for (UInt i = 0; i<cc->getLayerCount(); ++i)
    {
    	const LayerData& layer = cc->getLayer(i);
    	//add item
    	item = new QListWidgetItem( layer_manager_ );
			name = layer.name.toQString();
			if (layer.flipped)
			{
				name += " (flipped)";
			}
			item->setText(name);
    	if (layer.visible)
    	{
    		item->setCheckState(Qt::Checked);
    	}
    	else
    	{
    		item->setCheckState(Qt::Unchecked);
    	}
    	if (layer.modified)
    	{
    		item->setText(item->text() + '*');
    	}
    	//highlight active item
    	if (i == cc->activeLayerIndex())
    	{
				layer_manager_->setCurrentItem(item);
    	}
		}
		layer_manager_->blockSignals(false);
  }
  
  void TOPPViewBase::updateSpectrumBar()
  {		
  	SpectrumCanvas* cc = activeCanvas_();
  	int layer_row = layer_manager_->currentRow();
  	if (layer_row == -1 || cc == 0)
  	{
  		return;
  	}
  	
  	spectrum_selection_->clear();
  	spectrum_selection_->blockSignals(true);
  	const LayerData& cl = cc->getCurrentLayer();
  	QTreeWidgetItem* item = 0;

  	if(cl.type == LayerData::DT_PEAK)
  	{
  		std::vector<QTreeWidgetItem*> parent_stack;
  		parent_stack.push_back(0);
  		bool fail = false;
  		
			for (UInt i = 0; i < cl.peaks.size(); ++i)
			{
				if (i > 0)
				{
					if (cl.peaks[i].getMSLevel() == cl.peaks[i-1].getMSLevel() + 1)
					{
						item = new QTreeWidgetItem(parent_stack.back());
						parent_stack.resize(parent_stack.size()+1);
					}
					else if (cl.peaks[i].getMSLevel() == cl.peaks[i-1].getMSLevel())
					{
						if (parent_stack.size() == 1)
						{
							item = new QTreeWidgetItem((QTreeWidget*)0);
						}
						else
						{
							item = new QTreeWidgetItem(*(parent_stack.end()-2));
						}
					}
					else if (cl.peaks[i].getMSLevel() < cl.peaks[i-1].getMSLevel())
					{
						int level_diff = cl.peaks[i-1].getMSLevel() - cl.peaks[i].getMSLevel();
						int parent_index = 0;
						QTreeWidgetItem* parent = 0;
						if (parent_stack.size() - level_diff >= 2)
						{
							parent_index = parent_stack.size() - level_diff - 1;
							parent = parent_stack[parent_index];

							item = new QTreeWidgetItem(parent, parent_stack[parent_index+1]);
						}
						else
						{
							item = new QTreeWidgetItem((QTreeWidget*)0);
						}
						parent_stack.resize(parent_index+1);
					}
					else
					{
						std::cerr << "Cannot build treelike view for spectrum browser, generating flat list instead." << std::endl;
						fail = true;
						break;
					}
				}
				else
				{
					item = new QTreeWidgetItem((QTreeWidget*)0);
				}
				
				parent_stack.back() = item;
				if (parent_stack.size() == 1)
				{
					spectrum_selection_->addTopLevelItem(item);
				}
				
				item->setText(0, QString("MS") + QString::number(cl.peaks[i].getMSLevel()));
				item->setText(1, QString::number(cl.peaks[i].getRT()));
				item->setText(2, QString::number(cl.peaks[i].getPrecursorPeak().getPosition()[0]));
				item->setText(3, QString::number(i));
			}
			if (fail)
			{
				// generate flat list instead
				spectrum_selection_->clear();
				for (UInt i = 0; i < cl.peaks.size(); ++i)
				{
					item = new QTreeWidgetItem((QTreeWidget*)0);
					item->setText(0, QString("MS") + QString::number(cl.peaks[i].getMSLevel()));
					item->setText(1, QString::number(cl.peaks[i].getRT()));
					item->setText(2, QString::number(cl.peaks[i].getPrecursorPeak().getPosition()[0]));
					item->setText(3, QString::number(i));
					spectrum_selection_->addTopLevelItem(item);
				}
			}
  	}
  	else
  	{
  		item = new QTreeWidgetItem((QTreeWidget*)0);
  		item->setText(0, QString("Feature map"));
  		item->setText(1, QString("-"));
  		item->setText(2, QString("-"));
  		item->setText(3, QString::number(0));
			item->setFlags(!Qt::ItemIsEnabled);
			spectrum_selection_->addTopLevelItem(item);
			return; // leave signals blocked
  	}
  	
  	if (cl.peaks.size() == 1)
  	{
  		item->setFlags(!Qt::ItemIsEnabled);
  		return; // leave signals blocked
  	}
  	
  	spectrum_selection_->blockSignals(false);
  }

	void TOPPViewBase::layerSelectionChange(int i)
	{
		if (i!=-1)
		{
			activeCanvas_()->activateLayer(i);
			updateFilterBar();
			updateSpectrumBar();
		}
	}
	
	void TOPPViewBase::spectrumSelectionChange(QTreeWidgetItem* item, int /*column*/)
	{
		SpectrumCanvas* cc = activeCanvas_();
		const LayerData& cl = cc->getCurrentLayer();
		
		int index = item->text(3).toInt();
		
		FeatureMapType f_dummy;
		ConsensusMapType c_dummy;
		ExperimentType exp;
		exp.resize(1);
		exp[0] = cl.peaks[index];
		addData_(f_dummy, c_dummy, exp, false, false, true, cl.filename, cl.name + " (" + QString::number(cl.peaks[index].getRT()) + ")");
			
		//updateSpectrumBar();
	}

	void TOPPViewBase::layerContextMenu(const QPoint & pos)
	{
		QListWidgetItem* item = layer_manager_->itemAt(pos);
		if (item)
		{
			QAction* new_action = 0;
			int layer = layer_manager_->row(item);
			QMenu* context_menu = new QMenu(layer_manager_);
			context_menu->addAction("Rename");
			context_menu->addAction("Delete");
			
			if (activeCanvas_()->getLayer(layer).flipped)
			{
				new_action = context_menu->addAction("Flip upwards (1D)");
			}
			else
			{
				new_action = context_menu->addAction("Flip downwards (1D)");
			}
			if (!active1DWindow_())
			{
				new_action->setEnabled(false);
			}
			QAction* selected = context_menu->exec(layer_manager_->mapToGlobal(pos));
			//delete layer
			if (selected!=0 && selected->text()=="Delete")
			{
				activeCanvas_()->removeLayer(layer);
			}
			//rename layer
			else if (selected!=0 && selected->text()=="Rename")
			{
				QString name = QInputDialog::getText(this,"Rename layer","Name:");
				if (name!="")
				{
					activeCanvas_()->setLayerName(layer, name);
				}
			}
			// flip layer up/downwards
			else if (selected != 0 && selected->text() == "Flip downwards (1D)")
			{
				activeCanvas_()->getLayer(layer).flipped = true;
				active1DWindow_()->canvas()->setMirrorModeActive(true);
			}
			else if (selected != 0 && selected->text() == "Flip upwards (1D)")
			{
				activeCanvas_()->getLayer(layer).flipped = false;
				bool b = active1DWindow_()->canvas()->flippedLayersExist();
				active1DWindow_()->canvas()->setMirrorModeActive(b);
			}
			
			//Update tab bar and window title
			if (activeCanvas_()->getLayerCount()!=0)
			{
				tab_bar_->setTabText(tab_bar_->currentIndex(), activeCanvas_()->getLayer(0).name.toQString());
				activeWindow_()->setWindowTitle(activeCanvas_()->getLayer(0).name.toQString());
			}
			else
			{
				tab_bar_->setTabText(tab_bar_->currentIndex(),"empty");
				activeWindow_()->setWindowTitle("empty");
			}
					
			//Update filter bar, spectrum bar and layer bar
			updateLayerBar();
			updateSpectrumBar();
			updateFilterBar();
			updateMenu();
			
			delete (context_menu);
		}
	}

	void TOPPViewBase::logContextMenu(const QPoint & pos)
	{
		QMenu* context_menu = new QMenu(log_);
		context_menu->addAction("Clear");

		QAction* selected = context_menu->exec(log_->mapToGlobal(pos));
		
		//clear text
		if (selected!=0 && selected->text()=="Clear")
		{
			log_->clear();
		}
		delete (context_menu);
	}

	void TOPPViewBase::filterContextMenu(const QPoint & pos)
	{
		//do nothing if no window is open
		if (activeCanvas_()==0) return;
		
		QMenu* context_menu = new QMenu(filters_);			

		//warn if the current layer is not visible
		String layer_name = String("Layer: ") + activeCanvas_()->getCurrentLayer().name;
		if (!activeCanvas_()->getCurrentLayer().visible)
		{
			layer_name += " (invisible)";
		}
		context_menu->addAction(layer_name.toQString())->setEnabled(false);
		context_menu->addSeparator();

		//add actions
		QListWidgetItem* item = filters_->itemAt(pos);
		if (item)
		{
			context_menu->addAction("Edit");
			context_menu->addAction("Delete");
		}
		else
		{
			context_menu->addAction("Add filter");
		}
		//results
		QAction* selected = context_menu->exec(filters_->mapToGlobal(pos));
		if (selected!=0)
		{
			if(selected->text()=="Delete")
			{
				DataFilters filters = activeCanvas_()->getCurrentLayer().filters;
				filters.remove(filters_->row(item));
				activeCanvas_()->setFilters(filters);
				updateFilterBar();
			}
			else if (selected->text()=="Edit")
			{
				filterEdit(item);
			}
			else if (selected->text()=="Add filter")
			{
				DataFilters filters = activeCanvas_()->getCurrentLayer().filters;
				DataFilters::DataFilter filter;
				DataFilterDialog dlg(filter, this);
				if (dlg.exec())
				{
					filters.add(filter);
					activeCanvas_()->setFilters(filters);
					updateFilterBar();
				}
			}
		}	
		delete (context_menu);
	}

	void TOPPViewBase::filterEdit(QListWidgetItem* item)
	{
		DataFilters filters = activeCanvas_()->getCurrentLayer().filters;
		DataFilters::DataFilter filter = filters[filters_->row(item)];
		DataFilterDialog dlg(filter, this);
		if (dlg.exec())
		{
			filters.replace(filters_->row(item),filter);
			activeCanvas_()->setFilters(filters);
			updateFilterBar();
		}
	}

	void TOPPViewBase::layerEdit(QListWidgetItem* /*item*/)
	{
		activeCanvas_()->showCurrentLayerPreferences();
	}

  void TOPPViewBase::updateFilterBar()
  {
  	//update filters
  	filters_->clear();
		
		SpectrumCanvas* canvas = activeCanvas_();
		if (canvas==0) return;
		if (canvas->getLayerCount()==0) return;
		
		const DataFilters& filters = activeCanvas_()->getCurrentLayer().filters;
		for (UInt i=0; i<filters.size(); ++i)
		{
			QListWidgetItem* item = new QListWidgetItem(filters_);
			item->setText(filters[i].toString().toQString());
		}
  
  	//update check box
  	filters_check_box_->setChecked(activeCanvas_()->getCurrentLayer().filters.isActive());
  }

	void TOPPViewBase::layerFilterVisibilityChange(bool on)
	{
		if (activeCanvas_())
		{
			activeCanvas_()->changeLayerFilterState(activeCanvas_()->activeLayerIndex(),on);
		}
	}

	void TOPPViewBase::layerVisibilityChange(QListWidgetItem* item)
	{
		int layer;
		bool visible;
		layer = layer_manager_->row(item);
		visible = activeCanvas_()->getLayer(layer).visible;
		
		if (item->checkState()==Qt::Unchecked && visible)
		{
			activeCanvas_()->changeVisibility(layer, false);
		}
		else if (item->checkState()==Qt::Checked && !visible)
		{
			activeCanvas_()->changeVisibility(layer, true);
		}
	}

  void TOPPViewBase::updateTabBar(QWidget* w)
  {
  	if (w)
  	{
  		Int window_id = qobject_cast<SpectrumWidget*>(w)->window_id;
  		tab_bar_->setCurrentId(window_id);
  	}
  }

  void TOPPViewBase::tileVertical()
  {
    // primitive horizontal tiling
    QWidgetList windows = ws_->windowList();
    if ( !windows.count() ) return;

    if (active1DWindow_()) active1DWindow_()->showNormal();
    if (active2DWindow_()) active2DWindow_()->showNormal();

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
    if ( !windows.count() ) return;

    if (active1DWindow_()) active1DWindow_()->showNormal();
    if (active2DWindow_()) active2DWindow_()->showNormal();

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
    connect(sw->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(updateToolBar()));
    connect(sw->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(updateSpectrumBar()));
    connect(sw->canvas(),SIGNAL(layerModficationChange(UInt,bool)),this,SLOT(updateLayerBar()));
    connect(sw,SIGNAL(sendStatusMessage(std::string,OpenMS::UInt)),this,SLOT(showStatusMessage(std::string,OpenMS::UInt)));
    connect(sw,SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
  
  	Spectrum2DWidget* sw2 = qobject_cast<Spectrum2DWidget*>(sw);
  	if (sw2 != 0)
  	{
  		connect(sw2->getHorizontalProjection(),SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
  		connect(sw2->getVerticalProjection(),SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
  		connect(sw2,SIGNAL(showSpectrumAs1D(int)),this,SLOT(showSpectrumAs1D(int)));
  	}
  	
	  sw->setWindowTitle(caption.toQString());
		
		//add tab with id  
  	static int window_counter = 4711;
  	sw->window_id = window_counter++;
			
    tab_bar_->addTab(caption.toQString(), sw->window_id);
    
    //connect slots and sigals for removing the widget from the bar, when it is closed
    //- through the menu entry
    //- through the tab bar
    //- thourgh the MDI close button
    connect(sw,SIGNAL(aboutToBeDestroyed(int)),tab_bar_,SLOT(removeId(int)));
    
    tab_bar_->setCurrentId(sw->window_id);

		//show first window maximized (only visible windows are in the list)
		if (ws_->windowList().count()==0)
		{
			sw->showMaximized();
		}
		else
		{
			sw->show();
		}
  }

  void TOPPViewBase::gotoDialog()
  {
  	activeWindow_()->showGoToDialog();
  }

  SpectrumWidget*  TOPPViewBase::activeWindow_() const
  {
  	if (!ws_->activeWindow()) return 0;
    return qobject_cast<SpectrumWidget*>(ws_->activeWindow());
  }
  
  SpectrumCanvas*  TOPPViewBase::activeCanvas_() const
  {
    SpectrumWidget* sw = qobject_cast<SpectrumWidget*>(ws_->activeWindow());
    if (sw == 0)
    {
    	return 0;
    }
    return sw->canvas();
  }

  Spectrum1DWidget* TOPPViewBase::active1DWindow_() const
  {
		Spectrum1DWidget* w = qobject_cast<Spectrum1DWidget*>(activeWindow_());
		if (!w) return 0;
		return w;
  }
  
  Spectrum2DWidget* TOPPViewBase::active2DWindow_() const
  {
		Spectrum2DWidget* w = qobject_cast<Spectrum2DWidget*>(activeWindow_());
		if (!w) return 0;
		return w;
  }

  Spectrum3DWidget* TOPPViewBase::active3DWindow_() const
  {
		Spectrum3DWidget* w = qobject_cast<Spectrum3DWidget*>(activeWindow_());
		if (!w) return 0;
		return w;
  }

  void TOPPViewBase::loadPreferences(String filename)
  {
    //compose default ini file path
    String default_ini_file = String(QDir::homePath()) + "/.TOPPView.ini";

    if (filename=="")
    {
      filename = default_ini_file;
    }

    //load preferences, if file exists
    if (File::exists(filename))
    {
    	bool error = false;
    	Param tmp;
    	tmp.load(filename);
    	//apply preferences if they are of the current TOPPView version
    	if(tmp.exists("preferences:version") && tmp.getValue("preferences:version").toString()==VersionInfo::getVersion())
    	{
    		try
    		{
      		setParameters(tmp);
    		}
    		catch (Exception::InvalidParameter& e)
    		{
    			error = true;
    		}
    	}
    	else
    	{
				error = true;
    	}
			//set parameters to defaults when something is fishy with the parameters file
			if (error)
			{
  			//reset parameters
  			setParameters(Param());

				cerr << "The TOPPView preferences files '" << filename << "' was ignored. It is no longer compatible with this TOPPView version and will be replaced." << endl;
			}
    }
    else if (filename != default_ini_file)
    {
    	cerr << "Unable to load INI File: '" << filename << "'" << endl;
    }
    param_.setValue("PreferencesFile" , filename);

    //set the recent files
    Param p = param_.copy("preferences:RecentFiles");
    if (p.size()!=0)
    {
      for (Param::ParamIterator it=p.begin() ; it!=p.end() ; ++it)
      {
      	QString filename = it->value.toQString();
      	if (File::exists(filename)) recent_files_.append(filename);
      }
    }
		
    updateRecentMenu_();
  }

  void TOPPViewBase::savePreferences()
  {
    // replace recent files
    param_.removeAll("preferences:RecentFiles");

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
      addDataFile(action->text(),true,true);
		}
	}

  QStringList TOPPViewBase::getFileList_()
  {
		String filter_all = "readable files (*.dta *.dta2d";
		String filter_single = "dta files (*.dta);;dta2d files (*.dta2d)";
#ifdef ANDIMS_DEF
		filter_all +=" *.cdf";
		filter_single += ";;ANDI/MS files (*.cdf)";
#endif
		filter_all += " *.mzML *.mzXML *.mzData *.featureXML *.consensusXML);;" ;
		filter_single +=";;mzML files (*.mzML);;mzXML files (*.mzXML);;mzData files (*.mzData);;feature map (*.featureXML);;consensus feature map (*.consensusXML);;XML files (*.xml);;all files (*)";
	
	 	return QFileDialog::getOpenFileNames(this, "Open file(s)", param_.getValue("preferences:default_path").toQString(), (filter_all+ filter_single).toQString());
  }

  void TOPPViewBase::openFileDialog()
  {
	 	QStringList files = getFileList_();
		for(QStringList::iterator it=files.begin();it!=files.end();it++)
		{
			addDataFile(*it,true,true);
		}
  }

	void TOPPViewBase::connectToDB_(DBConnection& db)
	{
		//get the password if unset
		if (!param_.exists("DBPassword"))
		{
			stringstream ss;
			ss << "Enter password for user '" << (String)param_.getValue("preferences:db:login") << "' at '"<< (String)param_.getValue("preferences:db:host")<<":"<<(String)param_.getValue("preferences:db:port")<<"' : ";
			bool ok;
			QString text = QInputDialog::getText(this, "TOPPView database password", ss.str().c_str(), QLineEdit::Password,QString::null, &ok);
			if ( ok )
			{
				param_.setValue("DBPassword",text);
			}
		}
		
		if (param_.exists("DBPassword"))
		{
			try
			{
				db.connect((String)param_.getValue("preferences:db:name"), (String)param_.getValue("preferences:db:login"),(String)param_.getValue("DBPassword"),(String)param_.getValue("preferences:db:host"),(UInt)param_.getValue("preferences:db:port"));
			}
			catch (DBConnection::InvalidQuery er)
			{
				param_.remove("DBPassword");
				showLogMessage_(LS_ERROR,"Unable to log in to the database server",String("Check the login data in the preferences!\nDatabase error message: ") + er.what());
			}
		}
	}
	
  void TOPPViewBase::openDatabaseDialog()
  {
		DBConnection db;
		connectToDB_(db);
		if (db.isConnected())
		{
			vector<UInt> result;	
			DBOpenDialog db_dialog(db,result,this);
			if (db_dialog.exec())
			{
				db.disconnect();
				for (vector<UInt>::iterator it = result.begin();it!=result.end();++it)
				{
 					addDataDB(*it,true);
				}
			}
		}
  }

	void TOPPViewBase::rerunTOPPTool()
	{
		//warn if hidden layer => wrong layer selected...
		const LayerData& layer = activeCanvas_()->getCurrentLayer();
		if (!layer.visible)
		{
  		showLogMessage_(LS_NOTICE,"The current layer is not visible","Have you selected the right layer for this action?");
		}
		
		//delete old input and output file
		File::remove(topp_.file_name + "_in");
		File::remove(topp_.file_name + "_out");

		//run the tool
		runTOPPTool_();
	}

	void TOPPViewBase::showTOPPDialog()
	{
		QAction* action = qobject_cast<QAction*>(sender());
		showTOPPDialog_(action->data().toBool());
	}

	void TOPPViewBase::showTOPPDialog_(bool visible)
	{
		//warn if hidden layer => wrong layer selected...
		const LayerData& layer = activeCanvas_()->getCurrentLayer();
		if (!layer.visible)
		{
  		showLogMessage_(LS_NOTICE,"The current layer is not visible","Have you selected the right layer for this action?");
		}
		
		//create and store unique file name prefix for files
		topp_.file_name = param_.getValue("preferences:tmp_file_path").toString() + "/TOPPView_" + File::getUniqueName();
		if (!File::writable(topp_.file_name+"_ini"))
		{
			showLogMessage_(LS_ERROR,"Cannot create temporary file",String("Cannot write to '")+topp_.file_name+"'_ini!");
			return;
		}
		ToolsDialog tools_dialog(this,topp_.file_name+"_ini",param_.getValue("preferences:default_path"),getCurrentLayer()->type);
		
		if(tools_dialog.exec()==QDialog::Accepted)
		{
			//Store tool name, input parameter and output parameter
			topp_.tool = tools_dialog.getTool();
			topp_.in = tools_dialog.getInput();
			topp_.out = tools_dialog.getOutput();
			topp_.visible = visible;
			//run the tool
			runTOPPTool_();
		}
	}

	void TOPPViewBase::runTOPPTool_()
	{
		const LayerData& layer = activeCanvas_()->getCurrentLayer();
		
		//test if files are writable
		if (!File::writable(topp_.file_name+"_in"))
		{
			showLogMessage_(LS_ERROR,"Cannot create temporary file",String("Cannot write to '")+topp_.file_name+"_in'!");
			return;
		}
		if (!File::writable(topp_.file_name+"_out"))
		{
			showLogMessage_(LS_ERROR,"Cannot create temporary file",String("Cannot write to '")+topp_.file_name+"'_out!");
			return;
		}
		
		//Store data
		topp_.layer_name = layer.name;
		topp_.window_id = activeWindow_()->window_id;
		if (layer.type==LayerData::DT_PEAK)
		{
			MzDataFile f;
			f.setLogType(ProgressLogger::GUI);
			if (topp_.visible)
			{
				ExperimentType exp;
				activeCanvas_()->getVisiblePeakData(exp);
				f.store(topp_.file_name+"_in",exp);
			}
			else
			{
				f.store(topp_.file_name+"_in",layer.peaks);
			}
		}
		else if (layer.type==LayerData::DT_FEATURE)
		{
			if (topp_.visible)
			{
				FeatureMapType map;
				activeCanvas_()->getVisibleFeatureData(map);
				FeatureXMLFile().store(topp_.file_name+"_in",map);
			}
			else
			{
				FeatureXMLFile().store(topp_.file_name+"_in",layer.features);
			}
		}
		else
		{
			if (topp_.visible)
			{
				ConsensusMapType map;
				activeCanvas_()->getVisibleConsensusData(map);
				ConsensusXMLFile().store(topp_.file_name+"_in",map);
			}
			else
			{
				ConsensusXMLFile().store(topp_.file_name+"_in",layer.consensus);
			}
		}
		
		//compose argument list
		QStringList args;
		args << "-ini"
				 << (topp_.file_name + "_ini").toQString()
				 << QString("-%1").arg(topp_.in.toQString())
				 << (topp_.file_name + "_in").toQString()
				 << "-no_progress";
		if (topp_.out!="")
		{
			args << QString("-%1").arg(topp_.out.toQString())
					 << (topp_.file_name+"_out").toQString();
		}
		
		//start log and show it
		showLogMessage_(LS_NOTICE,"Starting TOPP tool","");// tool + args.join(" "));
		
		//start process
		topp_.process = new QProcess();
		topp_.process->setProcessChannelMode(QProcess::MergedChannels);
		connect(topp_.process,SIGNAL(readyReadStandardOutput()),this,SLOT(updateProcessLog()));
		topp_.process->start(topp_.tool.toQString(),args);
		
		//connect the finished slot
		connect(topp_.process,SIGNAL(finished(int,QProcess::ExitStatus)),this,SLOT(finishTOPPToolExecution(int,QProcess::ExitStatus)));
		
		//start process
		topp_.process->waitForStarted();
		updateMenu();
	}

  void TOPPViewBase::finishTOPPToolExecution(int, QProcess::ExitStatus)
  {
  	//finish with new line
  	log_->append("");
  	
  	String tmp_dir = param_.getValue("preferences:tmp_file_path").toString();
  	
		if (topp_.process->exitStatus()==QProcess::CrashExit)
		{
			showLogMessage_(LS_ERROR,"Execution of TOPP tool not successful!",String("The tool crashed during execution. If you want to debug this crash, check the input files in '") + tmp_dir + "' or enable 'debug' mode in the TOPP ini file.");
		}
		else if(topp_.out!="")
		{
			if (!File::readable(topp_.file_name+"_out"))
			{
				showLogMessage_(LS_ERROR,"Cannot read TOPP output",String("Cannot read '")+topp_.file_name+"_out'!");
			}
			else
			{
				addDataFile(topp_.file_name+"_out",true,false, topp_.layer_name + " (" + topp_.tool + ")", topp_.window_id);
			}
		}
		
		//clean up
		delete topp_.process;
		topp_.process = 0;
		updateMenu();
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

  void TOPPViewBase::toggleProjections()
  {
    Spectrum2DWidget* w = active2DWindow_();
    if (w)
    {
    	w->toggleProjections();
    }
  }

	void TOPPViewBase::annotateWithID()
	{
		const LayerData& layer = activeCanvas_()->getCurrentLayer();
		//warn if hidden layer => wrong layer selected...
		if (!layer.visible)
		{
  		showLogMessage_(LS_NOTICE,"The current layer is not visible","Have you selected the right layer for this action?");
		}
				
		//load id data
		QString name = QFileDialog::getOpenFileName(this,"Select ProteinIdentification data",param_.getValue("preferences:default_path").toQString(),"identfication files (*.idXML);; all files (*.*)");
		if(name!="")
		{
			vector<PeptideIdentification> identifications; 
			vector<ProteinIdentification> protein_identifications; 
			IdXMLFile().load(name, protein_identifications, identifications);
			if (layer.type==LayerData::DT_PEAK)
			{
				IDMapper().annotate(const_cast<LayerData&>(layer).peaks, identifications, protein_identifications);
			}
			else if (layer.type==LayerData::DT_FEATURE)
			{
				IDMapper().annotate(const_cast<LayerData&>(layer).features,identifications,protein_identifications);
			}
			else
			{
				IDMapper().annotate(const_cast<LayerData&>(layer).consensus,identifications,protein_identifications);
			}
		}
	}

	void TOPPViewBase::showSpectrumGenerationDialog()
	{
		TheoreticalSpectrumGenerationDialog spec_gen_dialog;
		if (spec_gen_dialog.exec())
		{
			String seq_string(spec_gen_dialog.line_edit->text());
			if (seq_string == "")
			{
				QMessageBox::warning(this, "Error", "You must enter a peptide sequence!");
				return;
			}
			AASequence aa_sequence(seq_string);
			
			Int charge = spec_gen_dialog.spin_box->value();
			
			if (aa_sequence.isValid())
			{
				RichPeakSpectrum rich_spec;
				TheoreticalSpectrumGenerator generator;
				Param p;

				bool losses = (spec_gen_dialog.list_widget->item(7)->checkState() == Qt::Checked); // "Neutral losses"
				p.setValue("add_losses", losses, "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
				bool isotopes = (spec_gen_dialog.list_widget->item(8)->checkState() == Qt::Checked); // "Isotope clusters"
				p.setValue("add_isotopes", isotopes, "If set to 1 isotope peaks of the product ion peaks are added");
				p.setValue("a_intensity", spec_gen_dialog.a_intensity->value(), "Intensity of the a-ions");
				p.setValue("b_intensity", spec_gen_dialog.b_intensity->value(), "Intensity of the b-ions");
				p.setValue("c_intensity", spec_gen_dialog.c_intensity->value(), "Intensity of the c-ions");
				p.setValue("x_intensity", spec_gen_dialog.x_intensity->value(), "Intensity of the x-ions");
				p.setValue("y_intensity", spec_gen_dialog.y_intensity->value(), "Intensity of the y-ions");
				p.setValue("z_intensity", spec_gen_dialog.z_intensity->value(), "Intensity of the z-ions");
				DoubleReal rel_loss_int = (DoubleReal)(spec_gen_dialog.rel_loss_intensity->value()) / 100.0;
				p.setValue("relative_loss_intensity", rel_loss_int, "Intensity of loss ions, in relation to the intact ion intensity");
				generator.setParameters(p);

				if (spec_gen_dialog.list_widget->item(0)->checkState() == Qt::Checked) // "A-ions"
				{
					generator.addPeaks(rich_spec, aa_sequence, Residue::AIon, charge);
				}
				if (spec_gen_dialog.list_widget->item(1)->checkState() == Qt::Checked) // "B-ions"
				{
					generator.addPeaks(rich_spec, aa_sequence, Residue::BIon, charge);
				}
				if (spec_gen_dialog.list_widget->item(2)->checkState() == Qt::Checked) // "C-ions"
				{
					generator.addPeaks(rich_spec, aa_sequence, Residue::CIon, charge);
				}
				if (spec_gen_dialog.list_widget->item(3)->checkState() == Qt::Checked) // "X-ions"
				{
					generator.addPeaks(rich_spec, aa_sequence, Residue::XIon, charge);
				}
				if (spec_gen_dialog.list_widget->item(4)->checkState() == Qt::Checked) // "Y-ions"
				{
					generator.addPeaks(rich_spec, aa_sequence, Residue::YIon, charge); 
				}
				if (spec_gen_dialog.list_widget->item(5)->checkState() == Qt::Checked) // "Z-ions"
				{
					generator.addPeaks(rich_spec, aa_sequence, Residue::ZIon, charge); 
				}
				if (spec_gen_dialog.list_widget->item(6)->checkState() == Qt::Checked) // "Precursor"
				{
					generator.addPrecursorPeaks(rich_spec, aa_sequence, charge);
				}
				
				PeakSpectrum new_spec;
				for (RichPeakSpectrum::Iterator it = rich_spec.begin(); it != rich_spec.end(); ++it)
				{
					new_spec.push_back(static_cast<Peak1D>(*it));
				}
				
				PeakMap new_exp;
				new_exp.push_back(new_spec);
				
				FeatureMapType f_dummy;
				ConsensusMapType c_dummy;
				addData_(f_dummy, c_dummy, new_exp, false, false, true, "", seq_string + QString(" (theoretical)"));
	      
	      // ensure spectrum is drawn as sticks
	      draw_group_1d_->button(Spectrum1DCanvas::DM_PEAKS)->setChecked(true);
				setDrawMode1D(Spectrum1DCanvas::DM_PEAKS);
			}
			else
			{
				QMessageBox::warning(this, "Error", "The entered peptide sequence is invalid!");
			}
		}
	}
	
	void TOPPViewBase::showSpectrumAlignmentDialog()
	{
		SpectrumAlignmentDialog spec_align_dialog;
		if (spec_align_dialog.exec())
		{
// 			Spectrum1DMirrorWidget* active_1d_mirror_window = active1DMirrorWindow_();
// 			// only possible in mirror mode:
// 			if (active_1d_mirror_window)
// 			{
// 				SpectrumAlignment aligner;
// 				Param param;
// 				DoubleReal tolerance = spec_align_dialog.tolerance_spinbox->value();
// 				param.setValue("tolerance", tolerance, "Defines the absolut (in Da) or relative (in ppm) tolerance", false);
// 				String unit_is_ppm = spec_align_dialog.ppm->isChecked() ? "true" : "false";
// 				param.setValue("is_relative_tolerance", unit_is_ppm, "If true, the 'tolerance' is interpreted as ppm-value", false);
// 				aligner.setParameters(param);
// 				
// 				// TODO JJ
// 				
// 				const LayerData& current_layer_1 = active_1d_mirror_window->canvas()->getCurrentLayer();
// 				const LayerData& current_layer_2 = active_1d_mirror_window->flippedCanvas()->getCurrentLayer();
// 				const ExperimentType& map_1 = current_layer_1.peaks;
// 				const ExperimentType& map_2 = current_layer_2.peaks;
// 				const ExperimentType::SpectrumType& spectrum_1 = *(map_1.begin());
// 				const ExperimentType::SpectrumType& spectrum_2 = *(map_2.begin());
// 				std::vector<std::pair<UInt, UInt> > alignment;
// 
// 				aligner.getSpectrumAlignment(alignment, spectrum_1, spectrum_2);
// 				
// 				std::vector<std::pair<DoubleReal, DoubleReal > > alignment_lines;
// 				
// 				for (UInt i = 0; i < alignment.size(); ++i)
// 				{
// 					DoubleReal line_begin_mz = spectrum_1[alignment[i].first].getMZ();
// 					DoubleReal line_end_mz = spectrum_2[alignment[i].second].getMZ();
// 					alignment_lines.push_back(std::make_pair(line_begin_mz, line_end_mz));
// 				}
// 				active_1d_mirror_window->setAlignmentLines(alignment_lines);
// 				
// 				SpectrumAlignmentScore scorer;
// 				scorer.setParameters(param);
// 				double score = scorer(spectrum_1, spectrum_2);
// 				
// 				QMessageBox::information(this, "Alignment performed", QString("Aligned %1 pairs of peaks (Score: %2).").arg(alignment_lines.size()).arg(score));
// 			}
// 			else
// 			{
// 				QMessageBox::warning(this, "Not supported", "Here be some description.");
// 			}
		}
	}
	
	void TOPPViewBase::showCurrentPeaksAs3D()
	{
    const LayerData& layer = activeCanvas_()->getCurrentLayer();
  	if (layer.type==LayerData::DT_PEAK)
  	{
  		//open new 3D widget
  		Spectrum3DWidget* w = new Spectrum3DWidget(getSpectrumParameters_(3), ws_);
			
  		//copy data
  		ExperimentType exp;
			activeCanvas_()->getVisiblePeakData(exp);
  			
	    if (!w->canvas()->addLayer(exp))
	  	{
	  		return;
	  	}
			String caption = layer.name + " (3D)";
			w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);
      showAsWindow_(w,caption);
	    updateLayerBar();
	    updateSpectrumBar();
			updateFilterBar();
			updateMenu();	
		}
		else
		{
      showLogMessage_(LS_NOTICE,"Wrong layer type","You cannot open feature data in 3D mode.");
		}
	}

	void TOPPViewBase::showSpectrumAs1D(int index)
	{
		const LayerData& layer = activeCanvas_()->getCurrentLayer();
		//copy spectrum
		ExperimentType exp;
		exp.resize(1);
		exp[0] = layer.peaks[index];
		//open new 1D widget
		Spectrum1DWidget* w = new Spectrum1DWidget(getSpectrumParameters_(1), ws_);
    
    //add data
    if (!w->canvas()->addLayer(exp))
  	{
  		return;
  	}
    
		String caption = layer.name + " (RT: " + String(layer.peaks[index].getRT()) + ")";
		w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);
		
    showAsWindow_(w,caption);
    updateLayerBar();
    updateSpectrumBar();
		updateFilterBar();
		updateMenu();
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
													 "Lesser GNU Public License (LGPL)<BR>"
													 "<BR>"
													 "<BR>"
													 "<BR>"
													 "<BR>"
													 "<BR>"
													 "Any published work based on TOPP and OpenMS shall cite these papers:<BR>"
													 "Sturm et al., BMC Bioinformatics (2008), 9, 163<BR>"
													 "Kohlbacher et al., Bioinformatics (2007), 23:e191-e197<BR>"
													 ).arg(VersionInfo::getVersion().toQString());
		label = new QLabel(text,dlg);
		grid->addWidget(label,0,1,Qt::AlignTop | Qt::AlignLeft);
		
		//execute
		dlg->exec();
	}

	void TOPPViewBase::updateProcessLog()
	{
		//show log if there is output
		qobject_cast<QWidget *>(log_->parent())->show();
		
		//update log_
		log_->textCursor().insertText(topp_.process->readAllStandardOutput());
	}

	Param TOPPViewBase::getSpectrumParameters_(UInt dim)
	{
		Param out = param_.copy(String("preferences:") + dim + "d:",true);
		out.setValue("default_path",param_.getValue("preferences:default_path").toString());
		out.setValue("on_file_change",param_.getValue("preferences:on_file_change").toString());
		return out;
	}

  void TOPPViewBase::abortTOPPTool()
  {
  	if (topp_.process)
  	{
  		//block signals to avoid error message from finished() signal
  		topp_.process->blockSignals(true);
  		//kill and delete the process
  		topp_.process->terminate();
  		delete topp_.process;
  		topp_.process = 0;

			//finish log with new line
	  	log_->append("");

  		updateMenu();
  	}
  }
  
  void TOPPViewBase::updateMenu()
  {
  	//is there a canvas?
  	bool canvas_exists = false;
  	if (activeCanvas_()!=0)
  	{
  		canvas_exists = true;
  	}
  	//is there a layer?
  	bool layer_exists = false;
  	if (canvas_exists && activeCanvas_()->getLayerCount()!=0)
  	{
  		layer_exists = true;
  	}
		//is there a TOPP tool running
		bool topp_running = false;
		if (topp_.process!=0)
		{
			topp_running = true;
		}
  	bool mirror_mode = false;
  	if (active1DWindow_() && active1DWindow_()->canvas()->mirrorModeActive())
  	{
  		mirror_mode = true;
  	}
		QList<QAction*> actions = this->findChildren<QAction*>("");
		for (int i=0; i<actions.count(); ++i)
		{
			QString text = actions[i]->text();
			if (text=="&Close" || text=="Show/hide grid lines" || text=="Show/hide axis legends")
			{
				actions[i]->setEnabled(false);
				if (canvas_exists)
				{
					actions[i]->setEnabled(true);
				}
			}
			else if (text.left(15)=="Apply TOPP tool")
			{
				actions[i]->setEnabled(false);
				if (canvas_exists && layer_exists && !topp_running)
				{
					actions[i]->setEnabled(true);
				}
			}
			else if (text=="Abort running TOPP tool")
			{
				actions[i]->setEnabled(false);
				if (topp_running)
				{
					actions[i]->setEnabled(true);
				}
			}
			else if (text=="Rerun TOPP tool")
			{
				actions[i]->setEnabled(false);
				if (canvas_exists && layer_exists && !topp_running && topp_.tool!="")
				{
					actions[i]->setEnabled(true);
				}
			}
			else if (text=="&Go to" || text=="&Edit meta data" || text=="&Statistics" || text=="&Annotate with identifiction"  || text=="Save all data"  || text=="Save visible data"  || text=="Preferences")
			{
				actions[i]->setEnabled(false);
				if (canvas_exists && layer_exists)
				{
					actions[i]->setEnabled(true);
				}
			}
			else if (text=="Align spectra")
			{
				actions[i]->setEnabled(false);
				if (mirror_mode)
				{
					actions[i]->setEnabled(true);
				}
			}
			else if (text=="")
			{
				actions[i]->setEnabled(false);
				if (canvas_exists && layer_exists)
				{
					actions[i]->setEnabled(true);
				}
			}
		}
  }


  void TOPPViewBase::loadFiles(const StringList& list, QSplashScreen* splash_screen)
  {
		bool last_was_plus = false;
    for (StringList::const_iterator it=list.begin(); it!=list.end(); ++it)
    {
      if (*it=="+")
      {
      	last_was_plus = true;
      	continue;
    	}
    	else if (*it=="@bw")
    	{
    		if ( (active2DWindow_()!=0 || active3DWindow_()!=0) && activeCanvas_()!=0 )
    		{
    			Param tmp = activeCanvas_()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#ffffff;100,#000000");
    			activeCanvas_()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (*it=="@bg")
    	{
    		if ( (active2DWindow_()!=0 || active3DWindow_()!=0) && activeCanvas_()!=0 )
    		{
    			Param tmp = activeCanvas_()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#dddddd;100,#000000");
    			activeCanvas_()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (*it=="@b")
    	{
    		if ( (active2DWindow_()!=0 || active3DWindow_()!=0) && activeCanvas_()!=0 )
    		{
    			Param tmp = activeCanvas_()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#000000;100,#000000");
    			activeCanvas_()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (*it=="@r")
    	{
    		if ( (active2DWindow_()!=0 || active3DWindow_()!=0) && activeCanvas_()!=0 )
    		{
    			Param tmp = activeCanvas_()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#ff0000;100,#ff0000");
    			activeCanvas_()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (*it=="@g")
    	{
    		if ( (active2DWindow_()!=0 || active3DWindow_()!=0) && activeCanvas_()!=0 )
    		{
    			Param tmp = activeCanvas_()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#00ff00;100,#00ff00");
    			activeCanvas_()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (*it=="@m")
    	{
    		if ( (active2DWindow_()!=0 || active3DWindow_()!=0) && activeCanvas_()!=0 )
    		{
    			Param tmp = activeCanvas_()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#ff00ff;100,#ff00ff");
    			activeCanvas_()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (!last_was_plus || !activeWindow_())
    	{
    		splash_screen->showMessage((String("Loading file: ") + *it).toQString());
    		splash_screen->repaint();
    		QApplication::processEvents();
    		addDataFile(*it,false,true);
    	}
    	else 
    	{
    		splash_screen->showMessage((String("Loading file: ") + *it).toQString());
    		splash_screen->repaint();
    		QApplication::processEvents();
    		last_was_plus = false;
    		addDataFile(*it,false,true,"",activeWindow_()->window_id);
    	}
    }
  }

  void TOPPViewBase::showLogMessage_(TOPPViewBase::LogState state, const String& heading, const String& body)
  {
		//Compose current time string
		DateTime d;
		d.now();
		
		String state_string;
		switch(state)
		{
			case LS_NOTICE: state_string = "NOTICE"; break;
			case LS_WARNING: state_string = "WARNING"; break;
			case LS_ERROR: state_string = "ERROR"; break;
		}
		
		//update log
		log_->append("==============================================================================");
		log_->append((d.getTime() + " " + state_string + ": " + heading).toQString());
		log_->append(body.toQString());

  	//show log tool window
		qobject_cast<QWidget *>(log_->parent())->show();
  }

	void TOPPViewBase::saveLayerAll()
	{
		activeCanvas_()->saveCurrentLayer(false);
	}
	
	void TOPPViewBase::saveLayerVisible()
	{
		activeCanvas_()->saveCurrentLayer(true);		
	}
	
	void TOPPViewBase::toggleGridLines()
	{
		activeCanvas_()->showGridLines(!activeCanvas_()->gridLinesShown());		
	}
	
	void TOPPViewBase::toggleAxisLegends()
	{
		activeWindow_()->showLegend(!activeWindow_()->isLegendShown());		
	}
	
	void TOPPViewBase::showPreferences()
	{
		activeCanvas_()->showCurrentLayerPreferences();
	}
	
	void TOPPViewBase::metadataFileDialog()
	{
	 	QStringList files = getFileList_();
		FileHandler fh;
		fh.getOptions().setMetadataOnly(true);
		for(QStringList::iterator it=files.begin();it!=files.end();it++)
		{
			ExperimentType exp;
			try
			{
				fh.loadExperiment(*it,exp);
			}
			catch (Exception::BaseException& e)
			{
				QMessageBox::critical(this,"Error",(String("Error while reading data: ")+e.what()).c_str());
	      return;
			}
			MetaDataBrowser dlg(false, this);
			dlg.setWindowTitle("Meta data");			
			dlg.add(exp);
	 	 	dlg.exec();			
		}
	}

	void TOPPViewBase::metadataDatabaseDialog()
	{
		DBConnection con;
		connectToDB_(con);
		if (con.isConnected())
		{
			vector<UInt> ids;	
			DBOpenDialog db_dialog(con,ids,ws_);
			if (db_dialog.exec())
			{
				DBAdapter db(con);
				db.getOptions().setMetadataOnly(true);
				for (vector<UInt>::iterator it = ids.begin();it!=ids.end();++it)
				{
					ExperimentType exp;
					try
					{
						db.loadExperiment(*it, exp);
					}
					catch (Exception::BaseException& e)
					{
						QMessageBox::critical(this,"Error",(String("Error while reading data: ")+e.what()).c_str());
			      return;
					}
					MetaDataBrowser dlg(false, this);
					dlg.setWindowTitle("Meta data");			
					dlg.add(exp);
			 	 	dlg.exec();
				}
			}
		}
	}

	void TOPPViewBase::copyLayer(const QMimeData* /*data*/, int id)
	{
		//NOT USED RIGHT NOW, BUT KEEP THIS CODE (it was hard to find out how this is done)
		//decode data to get the row
		//QByteArray encoded_data = data->data(data->formats()[0]);
		//QDataStream stream(&encoded_data, QIODevice::ReadOnly);
		//int row, col;
		//stream >> row >> col;

  	//set wait cursor
  	setCursor(Qt::WaitCursor);		
		
		//only the selected row can be dragged => the source layer is the selected layer
		const LayerData& layer = activeCanvas_()->getCurrentLayer();
		
		//copy the feature and peak data
		FeatureMapType features = layer.features;
		ExperimentType peaks = layer.peaks;
		ConsensusMapType consensus = layer.consensus;
		
		//determine where to copy the data
		UInt new_id = 0;
		if (id!=-1) new_id = id;

		//determine if the data is 2D data
		bool is_2D = false;
		bool is_feature = false;
    if (layer.type==LayerData::DT_FEATURE)
    {
      is_2D = true;
      is_feature = true;
    }
    else if (layer.type==LayerData::DT_CONSENSUS)
    {
      is_2D = true;
      is_feature = true;
    }
    else
    {
    	UInt ms1_scans = 0;
    	for (UInt i=0; i<peaks.size();++i)
    	{
    		if (peaks[i].getMSLevel()==1) ++ms1_scans;
    		if (ms1_scans>1)
    		{
    			is_2D = true;
    			break;
    		}
    	}
    }
		
		//add the data
		addData_(features, consensus, peaks, is_feature, is_2D, false, layer.filename, layer.name, new_id);

		//reset cursor
  	setCursor(Qt::ArrowCursor);		
	}

	void TOPPViewBase::keyPressEvent(QKeyEvent* e)
	{
 		SpectrumCanvas* canvas = activeCanvas_();
    if (canvas == 0 || canvas->getLayerCount()==0)
    {
    	e->ignore();
      return;
    }
    
		//page up => go one layer up
		if (e->key()==Qt::Key_PageUp)
		{
			if (canvas->activeLayerIndex()!=0)
			{
				canvas->activateLayer(canvas->activeLayerIndex()-1);
				updateLayerBar();
				updateFilterBar();
				updateMenu();
				e->accept();
			}
		}
		//page down => go one layer down
		else if (e->key()==Qt::Key_PageDown)
		{
			if (canvas->activeLayerIndex()!=canvas->getLayerCount()-1)
			{
				canvas->activateLayer(canvas->activeLayerIndex()+1);
				updateLayerBar();
				updateFilterBar();
				updateMenu();
				e->accept();
			}
		}
		
		e->ignore();
	}

} //namespace OpenMS

