// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>

#include <OpenMS/CONCEPT/EnumHelpers.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DataSelectionTabs.h>
#include <OpenMS/VISUAL/DIALOGS/SpectrumAlignmentDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewPrefDialog.h>
#include <OpenMS/VISUAL/INTERFACES/IPeptideIds.h>
#include <OpenMS/VISUAL/LayerListView.h>
#include <OpenMS/VISUAL/LayerDataChrom.h>
#include <OpenMS/VISUAL/LayerDataConsensus.h>
#include <OpenMS/VISUAL/LayerDataFeature.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/LogWindow.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>
#include <OpenMS/VISUAL/Plot1DWidget.h>
#include <OpenMS/VISUAL/Plot2DCanvas.h>
#include <OpenMS/VISUAL/Plot2DWidget.h>
#include <OpenMS/VISUAL/Plot3DCanvas.h>
#include <OpenMS/VISUAL/Plot3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/Plot3DWidget.h>
#include <OpenMS/VISUAL/SpectraIDViewTab.h>
#include <OpenMS/VISUAL/SpectraTreeTab.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>

// Qt
#include <QCloseEvent>
#include <QtCore/QDir>
#include <QtCore/QSettings>
#include <QtCore/QUrl>
#include <QGuiApplication>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QSplashScreen>
#include <QtWidgets/QToolButton>

#include <cmath>
#include <utility>

using namespace std;

namespace OpenMS
{
  using namespace Internal;
  using namespace Math;

  const std::string user_section = "preferences:user:";

  /// supported types which can be opened with File-->Open
  const FileTypeList supported_types({ FileTypes::MZML, FileTypes::MZXML, FileTypes::MZDATA, FileTypes::SQMASS,
                                       FileTypes::FEATUREXML, FileTypes::CONSENSUSXML, FileTypes::IDXML,
                                       FileTypes::DTA, FileTypes::DTA2D, FileTypes::MGF, FileTypes::MS2,
                                       FileTypes::MSP, FileTypes::BZ2, FileTypes::GZ });

  TOPPViewBase::TOPPViewBase(TOOL_SCAN scan_mode, VERBOSITY verbosity, QWidget* parent) :
    QMainWindow(parent),
    DefaultParamHandler("TOPPViewBase"),
    scan_mode_(scan_mode),
    verbosity_(verbosity),
    ws_(this),
    tab_bar_(this),
    recent_files_(),
    menu_(this, &ws_, &recent_files_)
  {
    setWindowTitle("TOPPView");
    setWindowIcon(QIcon(":/TOPPView.png"));
    setMinimumSize(400, 400); // prevents errors caused by too small width, height values
    setAcceptDrops(true); // enable drag-and-drop

    // get geometry of first screen
    QRect screen_geometry = QGuiApplication::primaryScreen()->geometry();
    // center main window
    setGeometry(
      (int)(0.1 * screen_geometry.width()),
      (int)(0.1 * screen_geometry.height()),
      (int)(0.8 * screen_geometry.width()),
      (int)(0.8 * screen_geometry.height())
      );

    //################## Main Window #################
    // Create main workspace using a QVBoxLayout ordering the items vertically
    // (the tab bar and the main workspace). Uses a dummy central widget (to be able to
    // have a layout), then adds vertically the tab bar and workspace.
    QWidget* dummy_cw = new QWidget(this);
    setCentralWidget(dummy_cw);
    QVBoxLayout* box_layout = new QVBoxLayout(dummy_cw);

    // create empty tab bar and workspace which will hold the main visualization widgets (e.g. spectrawidgets...)
    tab_bar_.setWhatsThis("Tab bar<BR><BR>Close tabs through the context menu or by double-clicking them.<BR>The tab bar accepts drag-and-drop from the layer bar.");
    tab_bar_.addTab("dummy", 4710);
    tab_bar_.setMinimumSize(tab_bar_.sizeHint());
    tab_bar_.removeId(4710);
    connect(&tab_bar_, &EnhancedTabBar::currentIdChanged, this, &TOPPViewBase::showWindow);
    connect(&tab_bar_, &EnhancedTabBar::closeRequested, this, &TOPPViewBase::closeByTab);
    connect(&tab_bar_, &EnhancedTabBar::dropOnWidget, [this](const QMimeData* data, QWidget* source){ copyLayer(data, source); });
    connect(&tab_bar_, &EnhancedTabBar::dropOnTab, this, &TOPPViewBase::copyLayer);
    box_layout->addWidget(&tab_bar_);

    // Trigger updates only when the active subWindow changes and update it
    connect(&ws_, &EnhancedWorkspace::subWindowActivated, [this](QMdiSubWindow* window) {
      if (window && lastActiveSubwindow_ != window) /* 0 upon terminate */ updateBarsAndMenus();
      lastActiveSubwindow_ = window;
    });
    connect(&ws_, &EnhancedWorkspace::dropReceived, this, &TOPPViewBase::copyLayer);
    box_layout->addWidget(&ws_);

    //################## STATUS #################
    // create status bar
    message_label_ = new QLabel(statusBar());
    statusBar()->addWidget(message_label_, 1);

    x_label_ = new QLabel("RT: 12345678", statusBar());
    x_label_->setMinimumSize(x_label_->sizeHint());
    x_label_->setText("");
    statusBar()->addPermanentWidget(x_label_, 0);
    y_label_ = new QLabel("m/z: 123456780912", statusBar());
    y_label_->setMinimumSize(y_label_->sizeHint());
    y_label_->setText("");
    statusBar()->addPermanentWidget(y_label_, 0);

    //################## TOOLBARS #################
    // create toolbars and connect signals
    QToolButton* b;

    //--Basic tool bar for all views--
    tool_bar_ = addToolBar("Basic tool bar");
    tool_bar_->setObjectName("tool_bar");

    // intensity modes
    intensity_button_group_ = new QButtonGroup(tool_bar_);
    intensity_button_group_->setExclusive(true);

    b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/lin.png"));
    b->setToolTip("Intensity: Normal");
    b->setShortcut(Qt::Key_N);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Normal<BR><BR>Intensity is displayed unmodified.<BR>(Hotkey: N)");
    intensity_button_group_->addButton(b, PlotCanvas::IM_NONE);
    tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/percentage.png"));
    b->setToolTip("Intensity: Percentage");
    b->setShortcut(Qt::Key_P);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Percentage<BR><BR>Intensity is displayed as a percentage of the layer"
                    " maximum intensity. If only one layer is displayed this mode behaves like the"
                    " normal mode. If more than one layer is displayed intensities are aligned."
                    "<BR>(Hotkey: P)");
    intensity_button_group_->addButton(b, PlotCanvas::IM_PERCENTAGE);
    tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/snap.png"));
    b->setToolTip("Intensity: Snap to maximum displayed intensity");
    b->setShortcut(Qt::Key_S);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Snap to maximum displayed intensity<BR><BR> In this mode the"
                    " color gradient is adapted to the maximum currently displayed intensity."
                    "<BR>(Hotkey: S)");
    intensity_button_group_->addButton(b, PlotCanvas::IM_SNAP);
    tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/log.png"));
    b->setToolTip("Intensity: Use log scaling for colors");
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Logarithmic scaling of intensities for color calculation");
    intensity_button_group_->addButton(b, PlotCanvas::IM_LOG);
    tool_bar_->addWidget(b);

    connect(intensity_button_group_, &QButtonGroup::idClicked, this, &TOPPViewBase::setIntensityMode);
    tool_bar_->addSeparator();

    // common buttons
    QAction* reset_zoom_button = tool_bar_->addAction(QIcon(":/reset_zoom.png"), "Reset Zoom", this, &TOPPViewBase::resetZoom);
    reset_zoom_button->setWhatsThis("Reset zoom: Zooms out as far as possible and resets the zoom history.<BR>(Hotkey: Backspace)");

    tool_bar_->show();

    //--1D toolbar--
    tool_bar_1d_ = addToolBar("1D tool bar");
    tool_bar_1d_->setObjectName("1d_tool_bar");

    // draw modes 1D
    draw_group_1d_ = new QButtonGroup(tool_bar_1d_);
    draw_group_1d_->setExclusive(true);

    b = new QToolButton(tool_bar_1d_);
    b->setIcon(QIcon(":/peaks.png"));
    b->setToolTip("Peak mode");
    b->setShortcut(Qt::Key_I);
    b->setCheckable(true);
    b->setWhatsThis("1D Draw mode: Peaks<BR><BR>Peaks are displayed as sticks.");
    draw_group_1d_->addButton(b, Plot1DCanvas::DM_PEAKS);
    tool_bar_1d_->addWidget(b);

    b = new QToolButton(tool_bar_1d_);
    b->setIcon(QIcon(":/lines.png"));
    b->setToolTip("Raw data mode");
    b->setShortcut(Qt::Key_R);
    b->setCheckable(true);
    b->setWhatsThis("1D Draw mode: Raw data<BR><BR>Peaks are displayed as a continuous line.");
    draw_group_1d_->addButton(b, Plot1DCanvas::DM_CONNECTEDLINES);
    tool_bar_1d_->addWidget(b);

    connect(draw_group_1d_, &QButtonGroup::idClicked, this, &TOPPViewBase::setDrawMode1D);
    tool_bar_->addSeparator();

    //--2D peak toolbar--
    tool_bar_2d_peak_ = addToolBar("2D peak tool bar");
    tool_bar_2d_peak_->setObjectName("2d_tool_bar");

    dm_precursors_2d_ = tool_bar_2d_peak_->addAction(QIcon(":/precursors.png"), "Show fragment scan precursors");
    dm_precursors_2d_->setCheckable(true);
    dm_precursors_2d_->setWhatsThis("2D peak draw mode: Precursors<BR><BR>fragment scan precursor peaks are marked.<BR>(Hotkey: 1)");
    dm_precursors_2d_->setShortcut(Qt::Key_1);

    connect(dm_precursors_2d_, &QAction::toggled, this, &TOPPViewBase::changeLayerFlag);

    projections_2d_ = tool_bar_2d_peak_->addAction(QIcon(":/projections.png"), "Show Projections", this, &TOPPViewBase::toggleProjections);
    projections_2d_->setCheckable(true);
    projections_2d_->setWhatsThis("Projections: Shows projections of peak data along RT and MZ axis.<BR>(Hotkey: 2)");
    projections_2d_->setShortcut(Qt::Key_2);

    //--2D feature toolbar--
    tool_bar_2d_feat_ = addToolBar("2D feature tool bar");
    tool_bar_2d_feat_->setObjectName("2d_feature_tool_bar");

    dm_hull_2d_ = tool_bar_2d_feat_->addAction(QIcon(":/convexhull.png"), "Show feature convex hull");
    dm_hull_2d_->setCheckable(true);
    dm_hull_2d_->setWhatsThis("2D feature draw mode: Convex hull<BR><BR>The convex hull of the feature is displayed.<BR>(Hotkey: 5)");
    dm_hull_2d_->setShortcut(Qt::Key_5);
    connect(dm_hull_2d_, &QAction::toggled, this, &TOPPViewBase::changeLayerFlag);

    dm_hulls_2d_ = tool_bar_2d_feat_->addAction(QIcon(":/convexhulls.png"), "Show feature convex hulls");
    dm_hulls_2d_->setCheckable(true);
    dm_hulls_2d_->setWhatsThis("2D feature draw mode: Convex hulls<BR><BR>The convex hulls of the feature are displayed: One for each mass trace.<BR>(Hotkey: 6)");
    dm_hulls_2d_->setShortcut(Qt::Key_6);
    connect(dm_hulls_2d_, &QAction::toggled, this, &TOPPViewBase::changeLayerFlag);

    // feature labels:
    dm_label_2d_ = new QToolButton(tool_bar_2d_feat_);
    dm_label_2d_->setPopupMode(QToolButton::MenuButtonPopup);
    QAction* action2 = new QAction(QIcon(":/labels.png"), "Show feature annotation", dm_label_2d_);
    action2->setCheckable(true);
    action2->setWhatsThis("2D feature draw mode: Labels<BR><BR>Display different kinds of annotation next to features.<BR>(Hotkey: 7)");
    action2->setShortcut(Qt::Key_7);
    dm_label_2d_->setDefaultAction(action2);
    tool_bar_2d_feat_->addWidget(dm_label_2d_);
    connect(dm_label_2d_, &QToolButton::triggered, this, &TOPPViewBase::changeLabel);
    // button menu
    group_label_2d_ = new QActionGroup(dm_label_2d_);
    QMenu* menu = new QMenu(dm_label_2d_);
    for (Size i = 0; i < LayerDataBase::SIZE_OF_LABEL_TYPE; ++i)
    {
      QAction* temp = group_label_2d_->addAction(
        QString(LayerDataBase::NamesOfLabelType[i].c_str()));
      temp->setCheckable(true);
      if (i == 0) temp->setChecked(true);
      menu->addAction(temp);
    }
    dm_label_2d_->setMenu(menu);

    // unassigned peptide identifications:
    dm_unassigned_2d_ = new QToolButton(tool_bar_2d_feat_);
    dm_unassigned_2d_->setPopupMode(QToolButton::MenuButtonPopup);
    QAction* action_unassigned = new QAction(QIcon(":/unassigned.png"), "Show unassigned peptide identifications", dm_unassigned_2d_);
    action_unassigned->setCheckable(true);
    action_unassigned->setWhatsThis("2D feature draw mode: Unassigned peptide identifications<BR><BR>Show unassigned peptide identifications by precursor m/z or by peptide mass.<BR>(Hotkey: 8)");
    action_unassigned->setShortcut(Qt::Key_8);
    dm_unassigned_2d_->setDefaultAction(action_unassigned);
    tool_bar_2d_feat_->addWidget(dm_unassigned_2d_);
    connect(dm_unassigned_2d_, &QToolButton::triggered, this, &TOPPViewBase::changeUnassigned);
    // button menu
    group_unassigned_2d_ = new QActionGroup(dm_unassigned_2d_);
    menu = new QMenu(dm_unassigned_2d_);
    StringList options = {"Don't show", "Show by precursor m/z", "Show by peptide mass", "Show label meta data"};
    for (const String& opt : options)
    {
      QAction* temp = group_unassigned_2d_->addAction(opt.toQString());
      temp->setCheckable(true);
      if (opt == options.front()) temp->setChecked(true);
      menu->addAction(temp);
    }
    dm_unassigned_2d_->setMenu(menu);

    //--2D consensus toolbar--
    tool_bar_2d_cons_ = addToolBar("2D peak tool bar");
    tool_bar_2d_cons_->setObjectName("2d_peak_tool_bar");

    dm_elements_2d_ = tool_bar_2d_cons_->addAction(QIcon(":/elements.png"), "Show consensus feature element positions");
    dm_elements_2d_->setCheckable(true);
    dm_elements_2d_->setWhatsThis("2D consensus feature draw mode: Elements<BR><BR>The individual elements that make up the  consensus feature are drawn.<BR>(Hotkey: 9)");
    dm_elements_2d_->setShortcut(Qt::Key_9);
    connect(dm_elements_2d_, &QAction::toggled, this, &TOPPViewBase::changeLayerFlag);

    //--2D identifications toolbar--
    tool_bar_2d_ident_ = addToolBar("2D identifications tool bar");
    tool_bar_2d_ident_->setObjectName("2d_ident_tool_bar");

    dm_ident_2d_ = tool_bar_2d_ident_->addAction(QIcon(":/peptidemz.png"), "Use theoretical peptide mass for m/z positions (default: precursor mass)");
    dm_ident_2d_->setCheckable(true);
    dm_ident_2d_->setWhatsThis("2D peptide identification draw mode: m/z source<BR><BR>Toggle between precursor mass (default) and theoretical peptide mass as source for the m/z positions of peptide identifications.<BR>(Hotkey: 5)");
    dm_ident_2d_->setShortcut(Qt::Key_5);
    connect(dm_ident_2d_, &QAction::toggled, this, &TOPPViewBase::changeLayerFlag);

    //################## Dock widgets #################
    // This creates the dock widgets: 
    // Layers, Views, Filters, and the Log dock widget on the bottom (by default hidden).

    // layer dock widget
    layer_dock_widget_ = new QDockWidget("Layers", this);
    layer_dock_widget_->setObjectName("layer_dock_widget");
    addDockWidget(Qt::RightDockWidgetArea, layer_dock_widget_);
    layers_view_ = new LayerListView(layer_dock_widget_);
    
    connect(layers_view_, &LayerListView::layerDataChanged, this, &TOPPViewBase::updateBarsAndMenus);
    layer_dock_widget_->setWidget(layers_view_);
    menu_.addWindowToggle(layer_dock_widget_->toggleViewAction());

    // Views dock widget
    views_dockwidget_ = new QDockWidget("Views", this);
    views_dockwidget_->setObjectName("views_dock_widget");
    addDockWidget(Qt::BottomDockWidgetArea, views_dockwidget_);
    selection_view_ = new DataSelectionTabs(views_dockwidget_, this);
    views_dockwidget_->setWidget(selection_view_);

    

    // add hide/show option to dock widget
    menu_.addWindowToggle(views_dockwidget_->toggleViewAction());

    // filter dock widget
    filter_dock_widget_ = new QDockWidget("Data filters", this);
    filter_dock_widget_->setObjectName("filter_dock_widget");
    addDockWidget(Qt::BottomDockWidgetArea, filter_dock_widget_);
    filter_list_ = new FilterList(filter_dock_widget_);
    connect(filter_list_, &FilterList::filterChanged, [&](const DataFilters& filter) {
      getActiveCanvas()->setFilters(filter);
    });
    filter_dock_widget_->setWidget(filter_list_);
    menu_.addWindowToggle(filter_dock_widget_->toggleViewAction());

    // log window
    QDockWidget* log_bar = new QDockWidget("Log", this);
    log_bar->setObjectName("log_bar");
    addDockWidget(Qt::BottomDockWidgetArea, log_bar);
    log_ = new LogWindow(log_bar);
    log_bar->setWidget(log_);
    menu_.addWindowToggle(log_bar->toggleViewAction());

    // tabify dock widgets so they don't fill up the whole space
    QMainWindow::tabifyDockWidget(filter_dock_widget_, log_bar);
    QMainWindow::tabifyDockWidget(log_bar, views_dockwidget_);

    //################## DEFAULTS #################
    initializeDefaultParameters_();

    // store defaults in param_
    defaultsToParam_();

    // load param file
    loadPreferences();

    // set current path
    current_path_ = param_.getValue(user_section + "default_path").toString();

    // set plugin search path, create it if it does not already exist
    if (verbosity_ == VERBOSITY::VERBOSE) 
    {
      tool_scanner_.setVerbose(1);
    }

    String plugin_path = String(param_.getValue(user_section + "plugins_path").toString());
    tool_scanner_.setPluginPath(plugin_path, true);

    // update the menu
    updateMenu();

    connect(&recent_files_, &RecentFilesMenu::recentFileClicked, this, &TOPPViewBase::openFile);

    // restore window positions
    QSettings settings("OpenMS", "TOPPView");
    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());

    //######################### File System Watcher ###########################################
    watcher_ = new FileWatcher(this);
    connect(watcher_, &FileWatcher::fileChanged, this, &TOPPViewBase::fileChanged_);
  }

  void TOPPViewBase::initializeDefaultParameters_()
  {
    // FIXME: these parameters are declared again in TOPPViewPrefDialog, incl. their allowed values
    //        There should be one place to do this. E.g. generate the GUI automatically from a Param (or simply use ParamEditor for the whole thing)

    // all parameters in preferences:user: can be edited by the user in the preferences dialog

    // general
    defaults_.setValue(user_section + "default_map_view", "2d", "Default visualization mode for maps.");
    defaults_.setValidStrings(user_section + "default_map_view", {"2d","3d"});
    defaults_.setValue(user_section + "default_path", ".", "Default path for loading and storing files.");
    defaults_.setValue(user_section + "default_path_current", "true", "If the current path is preferred over the default path.");
    defaults_.setValidStrings(user_section + "default_path_current", {"true","false"});
    defaults_.setValue(user_section + "plugins_path", File::getUserDirectory() + "OpenMS_Plugins", "Default path for loading Plugins");
    defaults_.setValue(user_section + "intensity_cutoff", "off", "Low intensity cutoff for maps.");
    defaults_.setValidStrings(user_section + "intensity_cutoff", {"on","off"});
    defaults_.setValue(user_section + "on_file_change", "ask", "What action to take, when a data file changes. Do nothing, update automatically or ask the user.");
    defaults_.setValidStrings(user_section + "on_file_change", {"none","ask","update automatically"});
    defaults_.setValue(user_section + "use_cached_ms2", "false", "If possible, only load MS1 spectra into memory and keep MS2 spectra on disk (using indexed mzML).");
    defaults_.setValidStrings(user_section + "use_cached_ms2", {"true","false"});
    defaults_.setValue(user_section + "use_cached_ms1", "false", "If possible, do not load MS1 spectra into memory spectra into memory and keep MS2 spectra on disk (using indexed mzML).");
    defaults_.setValidStrings(user_section + "use_cached_ms1", {"true","false"});

    // FIXME: getCanvasParameters() depends on the exact naming of the param sections below!
    // 1d view
    defaults_.insert(user_section + "1d:", Plot1DCanvas(Param()).getDefaults());
    defaults_.setSectionDescription(user_section + "1d", "Settings for single spectrum view.");
    // 2d view
    defaults_.insert(user_section + "2d:", Plot2DCanvas(Param()).getDefaults());
    defaults_.setSectionDescription(user_section + "2d", "Settings for 2D map view.");
    // 3d view
    defaults_.insert(user_section + "3d:", Plot3DCanvas(Param()).getDefaults());
    defaults_.setSectionDescription(user_section + "3d", "Settings for 3D map view.");
    // identification view
    defaults_.insert(user_section + "idview:", SpectraIDViewTab(Param()).getDefaults());
    defaults_.setSectionDescription(user_section + "idview", "Settings for identification view.");

    // non-editable parameters

    // not in Dialog (yet?)
    defaults_.setValue("preferences:topp_cleanup", "true", "If the temporary files for calling of TOPP tools should be removed after the call.");
    defaults_.setValidStrings("preferences:topp_cleanup", {"true", "false"});

    defaults_.setValue("preferences:version", "none", "OpenMS version, used to check if the TOPPView.ini is up-to-date");
    subsections_.emplace_back("preferences:RecentFiles");

    // store defaults in param_
    defaultsToParam_();
  }

  void TOPPViewBase::closeEvent(QCloseEvent* event)
  {
    ws_.closeAllSubWindows();
    QSettings settings("OpenMS", "TOPPView");
    settings.setValue("geometry", saveGeometry());
    settings.setValue("windowState", saveState());
    event->accept();
  }
  
  void TOPPViewBase::preferencesDialog()
  {
    Internal::TOPPViewPrefDialog dlg(this);
    dlg.setParam(param_.copy(user_section, true));

    // --------------------------------------------------------------------
    // Execute dialog and update parameter object with user modified values
    if (dlg.exec())
    {
      param_.remove(user_section);
      param_.insert(user_section, dlg.getParam());
      savePreferences();
    }
  }

  TOPPViewBase::LOAD_RESULT TOPPViewBase::addDataFile(const String& filename, bool show_options, bool add_to_recent, String caption, UInt window_id, Size spectrum_id)
  {
    String abs_filename = File::absolutePath(filename);

    // check if the file exists
    if (!File::exists(abs_filename))
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Open file error", String("The file '") + abs_filename + "' does not exist!");
      return LOAD_RESULT::FILE_NOT_FOUND;
    }

    // determine file type
    FileHandler fh;
    FileTypes::Type file_type = fh.getType(abs_filename);
    if (file_type == FileTypes::UNKNOWN)
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Open file error", String("Could not determine file type of '") + abs_filename + "'!");
      return LOAD_RESULT::FILETYPE_UNKNOWN;
    }

    // abort if file type unsupported
    if (!supported_types.contains(file_type))
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Open file error", String("The type '") + FileTypes::typeToName(file_type) + "' is not supported in TOPPView!");
      return LOAD_RESULT::FILETYPE_UNSUPPORTED;
    }

    // try to load data and determine if it's 1D or 2D data

    // create shared pointer to main data types
    FeatureMapType* feature_map = new FeatureMapType();
    FeatureMapSharedPtrType feature_map_sptr(feature_map);

    ExperimentSharedPtrType peak_map_sptr(new ExperimentType());

    ConsensusMapType* consensus_map = new ConsensusMapType();
    ConsensusMapSharedPtrType consensus_map_sptr(consensus_map);

    vector<PeptideIdentification> peptides;
    // not needed in data but for auto annotation
    vector<ProteinIdentification> proteins;
    String annotate_path;

    LayerDataBase::DataType data_type(LayerDataBase::DT_UNKNOWN);

    ODExperimentSharedPtrType on_disc_peaks(new OnDiscMSExperiment);

    // lock the GUI - no interaction possible when loading...
    GUIHelpers::GUILock glock(this);

    bool cache_ms2_on_disc = (param_.getValue(user_section + "use_cached_ms2") == "true");
    bool cache_ms1_on_disc = (param_.getValue(user_section + "use_cached_ms1") == "true");

    try
    {
      if (file_type == FileTypes::FEATUREXML)
      {
        FileHandler().loadFeatures(abs_filename, *feature_map, {FileTypes::FEATUREXML});
        data_type = LayerDataBase::DT_FEATURE;
      }
      else if (file_type == FileTypes::CONSENSUSXML)
      {
        FileHandler().loadConsensusFeatures(abs_filename, *consensus_map, {FileTypes::CONSENSUSXML});
        data_type = LayerDataBase::DT_CONSENSUS;
      }
      else if (file_type == FileTypes::IDXML || file_type == FileTypes::MZIDENTML)
      {
        FileHandler().loadIdentifications(abs_filename, proteins, peptides, {FileTypes::IDXML, FileTypes::MZIDENTML});
        if (peptides.empty())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No peptide identifications found");
        }
        // check if RT (and sequence) information is present:
        vector<PeptideIdentification> peptides_with_rt;
        for (vector<PeptideIdentification>::const_iterator it =
               peptides.begin(); it != peptides.end(); ++it)
        {
          if (!it->getHits().empty() && it->hasRT())
          {
            peptides_with_rt.push_back(*it);
          }
        }
        Size diff = peptides.size() - peptides_with_rt.size();
        if (diff)
        {
          String msg = String(diff) + " peptide identification(s) without"
                                      " sequence and/or retention time information were removed.\n" +
                       peptides_with_rt.size() + " peptide identification(s) remaining.";
          log_->appendNewHeader(LogWindow::LogState::WARNING, "While loading file:", msg);
        }
        if (peptides_with_rt.empty())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No peptide identifications with sufficient information remaining.");
        }
        peptides.swap(peptides_with_rt);

        if (!proteins.empty())
        {
          StringList paths;
          proteins[0].getPrimaryMSRunPath(paths);

          for (const String &path : paths)
          {
            if (File::exists(path) && fh.getType(path) == FileTypes::MZML)
            {
              annotate_path = path;
            }
          }
          // annotation could not be found in file reference
          if (annotate_path.empty())
          {
            // try to find file with same path & name but with mzML extension
            auto target = fh.swapExtension(abs_filename, FileTypes::Type::MZML);
            if (File::exists(target))
            {
              annotate_path = target;
            }
          }

          if (!annotate_path.empty())
          {
            // open dialog for annotation on load
            QMessageBox msg_box;
            auto spectra_file_name = File::basename(annotate_path);
            msg_box.setText("Spectra data for identification data was found.");
            msg_box.setInformativeText(String("Annotate spectra in " + spectra_file_name + "?").toQString());
            msg_box.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
            msg_box.setDefaultButton(QMessageBox::Yes);
            auto ret = msg_box.exec();
            if (ret == QMessageBox::No)
            { // no annotation performed
              annotate_path = "";
            }
          }
        }
        data_type = LayerDataBase::DT_IDENT;
      }
      else
      {
        bool parsing_success = false;
        if (file_type == FileTypes::MZML)
        {
          // Load index only and check success (is it indexed?)
          Internal::IndexedMzMLHandler indexed_mzml_file;
          indexed_mzml_file.openFile(filename);
          if ( indexed_mzml_file.getParsingSuccess() && cache_ms2_on_disc)
          {
            // If it has an index, now load index and meta data
            on_disc_peaks->openFile(filename, false);
            OPENMS_LOG_INFO << "INFO: will use cached MS2 spectra" << std::endl;
            if (cache_ms1_on_disc)
            {
              OPENMS_LOG_INFO << "log INFO: will use cached MS1 spectra" << std::endl;
            }
            parsing_success = true;

            // Caching strategy: peak_map_sptr will contain a MSSpectrum entry
            // for each actual spectrum on disk. However, initially these will
            // only be populated by the meta data (all data except the actual
            // raw data) which will allow us to read out RT, MS level etc.
            //
            // In a second step (see below), we populate some of these maps
            // with actual spectra including raw data (allowing us to only
            // populate MS1 spectra with actual data).

            peak_map_sptr = on_disc_peaks->getMetaData();

            for (Size k = 0; k < indexed_mzml_file.getNrSpectra() && !cache_ms1_on_disc; k++)
            {
              if ( peak_map_sptr->getSpectrum(k).getMSLevel() == 1)
              {
                peak_map_sptr->getSpectrum(k) = on_disc_peaks->getSpectrum(k);
              }
            }
            for (Size k = 0; k < indexed_mzml_file.getNrChromatograms() && !cache_ms2_on_disc; k++)
            {
              peak_map_sptr->getChromatogram(k) = on_disc_peaks->getChromatogram(k);
            }

            // Load at least one spectrum into memory (TOPPView assumes that at least one spectrum is in memory)
            if (cache_ms1_on_disc && peak_map_sptr->getNrSpectra() > 0) peak_map_sptr->getSpectrum(0) = on_disc_peaks->getSpectrum(0);
          }
        }

        // Load all data into memory if e.g. other file type than mzML
        if (!parsing_success)
        {
          fh.loadExperiment(abs_filename, *peak_map_sptr, {file_type}, ProgressLogger::GUI, true, true);
        }
        OPENMS_LOG_INFO << "INFO: done loading all " << std::endl;

        // a mzML file may contain both, chromatogram and peak data
        // -> this is handled in PlotCanvas::addPeakLayer FIXME: No it's not!
        if (peak_map_sptr->getNrSpectra() > 0 && peak_map_sptr->getNrChromatograms() > 0)
        {
          OPENMS_LOG_WARN << "Your input data contains chromatograms and spectra, falling back to display spectra only." << std::endl;
          data_type = LayerDataBase::DT_PEAK;
        }
        else if (peak_map_sptr->getNrChromatograms() > 0)
        {
          data_type = LayerDataBase::DT_CHROMATOGRAM;
        }
        else if (peak_map_sptr->getNrSpectra() > 0)
        {
          data_type = LayerDataBase::DT_PEAK;
        }
        else
        {
          throw Exception::FileEmpty(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MzML filed doesn't have either spectra or chromatograms.");
        }
      }
    }
    catch (Exception::BaseException& e)
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Error while loading file:", e.what());
      return LOAD_RESULT::LOAD_ERROR;
    }

    // sort for m/z and update ranges of newly loaded data
    peak_map_sptr->sortSpectra(true);
    peak_map_sptr->updateRanges(1);

    // try to add the data
    if (caption == "")
    {
      caption = FileHandler::stripExtension(File::basename(abs_filename));
    }
    else
    {
      abs_filename = "";
    }

    glock.unlock();

    if (!annotate_path.empty())
    { // this opens a new window with raw data + annotation; we want the actual idXML data on top
      auto load_res = addDataFile(annotate_path, false, false);
      if (load_res == LOAD_RESULT::OK)
      {
        auto l = getCurrentLayer();
        if (l)
        {
          bool success = l->annotate(peptides, proteins);
          if (success)
          {
            log_->appendNewHeader(LogWindow::LogState::NOTICE, "Done", "Annotation finished. Open identification view to see results!");
          }
          else
          {
            log_->appendNewHeader(LogWindow::LogState::NOTICE, "Error", "Annotation failed.");
          }
        }
        window_id = getActivePlotWidget()->getWindowId(); // add ids on top of raw data
      }
    }

    addData(feature_map_sptr, 
      consensus_map_sptr, 
      peptides, 
      peak_map_sptr, 
      on_disc_peaks, 
      data_type, 
      false,   // show as 1D
      show_options, 
      true,    // as new window
      abs_filename, 
      caption, 
      window_id, 
      spectrum_id);

    // add to recent file
    if (add_to_recent)
    {
      addRecentFile_(filename);
    }

    // watch file contents for changes
    watcher_->addFile(abs_filename);

    return LOAD_RESULT::OK;
  }

  void TOPPViewBase::addData(const FeatureMapSharedPtrType& feature_map,
                             const ConsensusMapSharedPtrType& consensus_map,
                             vector<PeptideIdentification>& peptides,
                             const ExperimentSharedPtrType& peak_map,
                             const ODExperimentSharedPtrType& on_disc_peak_map,
                             LayerDataBase::DataType data_type,
                             bool show_as_1d,
                             bool show_options,
                             bool as_new_window,
                             const String& filename,
                             const String& caption,
                             UInt window_id,
                             Size spectrum_id)
  {
    // initialize flags with defaults from the parameters
    bool maps_as_2d = (param_.getValue(user_section + "default_map_view") == "2d");
    bool maps_as_1d = false;
    bool use_intensity_cutoff = (param_.getValue(user_section + "intensity_cutoff") == "on");
    bool is_dia_data = false;

    // feature, consensus feature and identifications can be merged
    bool mergeable = ((data_type == LayerDataBase::DT_FEATURE) ||
                      (data_type == LayerDataBase::DT_CONSENSUS) ||
                      (data_type == LayerDataBase::DT_IDENT));

    // only one peak spectrum? disable 2D as default
    if (peak_map->size() == 1) { maps_as_2d = false; }

    // set the window where (new layer) data could be opened in
    // get EnhancedTabBarWidget with given id
    EnhancedTabBarWidgetInterface* tab_bar_target = ws_.getWidget(window_id);

    // cast to PlotWidget
    PlotWidget* target_window = dynamic_cast<PlotWidget*>(tab_bar_target);

    if (tab_bar_target == nullptr)
    {
      target_window = getActivePlotWidget();
    }
    else
    {
      as_new_window = false;
    }

    // create dialog no matter if it is shown or not. It is used to determine the flags.
    TOPPViewOpenDialog dialog(caption, as_new_window, maps_as_2d, use_intensity_cutoff, this);

    // disable opening in new window when there is no active window or feature/ID data is to be opened, but the current window is a 3D window
    if (target_window == nullptr || (mergeable && dynamic_cast<Plot3DWidget*>(target_window)))
    {
      dialog.disableLocation(true);
    }

    // for feature/consensus/identification maps
    if (mergeable) 
    {
      dialog.disableDimension(true); // disable 1d/2d/3d option
      dialog.disableCutoff(false);

      // enable merge layers if a feature layer is opened and there are already features layers to merge it to
      if (target_window)
      {
        PlotCanvas* open_canvas = target_window->canvas();
        std::map<Size, String> layers;
        for (Size i = 0; i < open_canvas->getLayerCount(); ++i)
        {
          if (data_type == open_canvas->getLayer(i).type)
          {
            layers[i] = open_canvas->getLayer(i).getName();
          }
        }
        dialog.setMergeLayers(layers); // adds a dropdown
      }
    }
    
    // show options if requested
    if (show_options && !dialog.exec())
    {
      return;
    }
    as_new_window = dialog.openAsNewWindow();
    maps_as_2d = dialog.viewMapAs2D();
    maps_as_1d = dialog.viewMapAs1D();
    if (show_as_1d)
    {
      maps_as_1d = true;
      maps_as_2d = false;
    }

    use_intensity_cutoff = dialog.isCutoffEnabled();
    is_dia_data = dialog.isDataDIA();
    Int merge_layer = dialog.getMergeLayer();

    // If we are dealing with DIA data, store this directly in the peak map
    // (ensures we will keep track of this flag from now on).
    if (is_dia_data)
    {
      peak_map->setMetaValue("is_dia_data", "true");
    }

    // determine the window to open the data in
    if (as_new_window) // new window
    {
      if (maps_as_1d) // 2d in 1d window
      {
        target_window = new Plot1DWidget(getCanvasParameters(1), DIM::Y, &ws_);
      }
      else if (maps_as_2d || mergeable) // 2d or features/IDs
      {
        target_window = new Plot2DWidget(getCanvasParameters(2), &ws_);
      }
      else // 3d
      {
        target_window = new Plot3DWidget(getCanvasParameters(3), &ws_);
      }
    }

    if (merge_layer == -1) // add new layer to the window
    {
      if (data_type == LayerDataBase::DT_FEATURE) // features
      {
        if (!target_window->canvas()->addLayer(feature_map, filename, caption))
        {
          return;
        }
      }
      else if (data_type == LayerDataBase::DT_CONSENSUS) // consensus features
      {
        if (! target_window->canvas()->addLayer(consensus_map, filename, caption))
          return;
      }
      else if (data_type == LayerDataBase::DT_IDENT)
      {
        if (! target_window->canvas()->addLayer(peptides, filename, caption))
          return;
      }
      else // peaks or chrom
      {
        if (data_type == LayerDataBase::DT_PEAK && ! target_window->canvas()->addPeakLayer(peak_map, on_disc_peak_map, filename, caption, use_intensity_cutoff))
        {
          return;
        }
        
        if (data_type == LayerDataBase::DT_CHROMATOGRAM &&
            !target_window->canvas()->addChromLayer(peak_map, on_disc_peak_map, filename, caption))
        {
          return;
        }
        
        Plot1DWidget* open_1d_window = dynamic_cast<Plot1DWidget*>(target_window);
        if (open_1d_window)
        {
          open_1d_window->canvas()->activateSpectrum(spectrum_id);
        }
      }
    }
    else // merge feature/ID data into feature layer
    {
      Plot2DCanvas* canvas = qobject_cast<Plot2DCanvas*>(target_window->canvas());
      if (data_type == LayerDataBase::DT_CONSENSUS)
      {
        canvas->mergeIntoLayer(merge_layer, consensus_map);
      }
      else if (data_type == LayerDataBase::DT_FEATURE)
      {
        canvas->mergeIntoLayer(merge_layer, feature_map);
      }
      else if (data_type == LayerDataBase::DT_IDENT)
      {
        canvas->mergeIntoLayer(merge_layer, peptides);
      }
      // combine layer names
      canvas->setLayerName(merge_layer, canvas->getLayerName(merge_layer) + " + " + caption);
    }

    if (as_new_window)
    {
      showPlotWidgetInWindow(target_window);
    }

    // enable spectra view tab (not required anymore since selection_view_.update() will decide automatically)
    //selection_view_->show(DataSelectionTabs::SPECTRA_IDX);
  }

  void TOPPViewBase::addRecentFile_(const String& filename)
  {
    recent_files_.add(filename);
  }

  void TOPPViewBase::openFile(const String& filename)
  {
    addDataFile(filename, true, true);
  }

  void TOPPViewBase::closeByTab(int id)
  {
    QWidget* w = dynamic_cast<QWidget*>(ws_.getWidget(id));
    if (w)
    {
      QMdiSubWindow* parent = qobject_cast<QMdiSubWindow*>(w->parentWidget());
      if (parent->close()) updateBarsAndMenus();
    }
  }

  void TOPPViewBase::showWindow(int id)
  {
    auto* sw = dynamic_cast<PlotWidget*>(ws_.getWidget(id));
    if (!sw) return;
    sw->setFocus(); // triggers layerActivated...
  }

  void TOPPViewBase::closeTab()
  {
    ws_.activeSubWindow()->close();
  }

  void TOPPViewBase::editMetadata()
  {
    PlotCanvas* canvas = getActiveCanvas();

    // warn if hidden layer => wrong layer selected...
    if (!canvas->getCurrentLayer().visible)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "The current layer is not visible", "Have you selected the right layer for this action?");
    }

    // show editable meta data dialog
    canvas->showMetaData(true);
  }

  void TOPPViewBase::layerStatistics() const
  {
    getActivePlotWidget()->showStatistics();
  }

  void TOPPViewBase::showStatusMessage(const string& msg, OpenMS::UInt time)
  {
    if (time == 0)
    {
      message_label_->setText(msg.c_str());
      statusBar()->update();
    }
    else
    {
      statusBar()->showMessage(msg.c_str(), time);
    }
  }

  void TOPPViewBase::showCursorStatus(const String& x, const String& y)
  {
    message_label_->setText("");
    x_label_->setText(x.toQString());
    y_label_->setText(y.toQString());
    statusBar()->update();
  }

  void TOPPViewBase::resetZoom() const
  {
    PlotWidget* w = getActivePlotWidget();
    if (w != nullptr)
    {
      w->canvas()->resetZoom();
    }
  }

  void TOPPViewBase::setIntensityMode(int index)
  {
    PlotWidget* w = getActivePlotWidget();
    if (w)
    {
      intensity_button_group_->button(index)->setChecked(true);
      w->setIntensityMode((OpenMS::PlotCanvas::IntensityModes)index);
    }
  }

  void TOPPViewBase::setDrawMode1D(int index) const
  {
    Plot1DWidget* w = getActive1DWidget();
    if (w)
    {
      w->canvas()->setDrawMode((OpenMS::Plot1DCanvas::DrawModes)index);
    }
  }

  void TOPPViewBase::changeLabel(QAction* action)
  {
    bool set = false;

    // label type is selected
    for (Size i = 0; i < LayerDataBase::SIZE_OF_LABEL_TYPE; ++i)
    {
      if (action->text().toStdString() == LayerDataBase::NamesOfLabelType[i])
      {
        getActive2DWidget()->canvas()->setLabel(LayerDataBase::LabelType(i));
        set = true;
      }
    }

    // button is simply pressed
    if (!set)
    {
      if (getActive2DWidget()->canvas()->getCurrentLayer().label == LayerDataBase::L_NONE)
      {
        getActive2DWidget()->canvas()->setLabel(LayerDataBase::L_INDEX);
        dm_label_2d_->menu()->actions()[1]->setChecked(true);
      }
      else
      {
        getActive2DWidget()->canvas()->setLabel(LayerDataBase::L_NONE);
        dm_label_2d_->menu()->actions()[0]->setChecked(true);
      }
    }

    updateToolBar();
  }

  void TOPPViewBase::changeUnassigned(QAction* action)
  {
    // mass reference is selected
    if (action->text().toStdString() == "Don't show")
    {
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::F_UNASSIGNED, false);
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::I_PEPTIDEMZ, false);
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::I_LABELS, false);
    }
    else if (action->text().toStdString() == "Show by precursor m/z")
    {
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::F_UNASSIGNED, true);
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::I_PEPTIDEMZ, false);
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::I_LABELS, false);
    }
    else if (action->text().toStdString() == "Show by peptide mass")
    {
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::F_UNASSIGNED, true);
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::I_PEPTIDEMZ, true);
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::I_LABELS, false);
    }
    else if (action->text().toStdString() == "Show label meta data")
    {
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::F_UNASSIGNED, true);
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::I_PEPTIDEMZ, false);
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::I_LABELS, true);
    }
    else // button is simply pressed
    {
      bool previous = getActive2DWidget()->canvas()->getLayerFlag(LayerDataBase::F_UNASSIGNED);
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::F_UNASSIGNED,
                                                  !previous);
      if (previous) // now: don't show
      {
        dm_unassigned_2d_->menu()->actions()[0]->setChecked(true);
      }
      else // now: show by precursor
      {
        dm_unassigned_2d_->menu()->actions()[1]->setChecked(true);
      }
      getActive2DWidget()->canvas()->setLayerFlag(LayerDataBase::I_PEPTIDEMZ, false);
    }

    updateToolBar();
  }

  void TOPPViewBase::changeLayerFlag(bool on)
  {
    QAction* action = qobject_cast<QAction*>(sender());
    if (Plot2DWidget* win = getActive2DWidget())
    {
      // peaks
      if (action == dm_precursors_2d_)
      {
        win->canvas()->setLayerFlag(LayerDataBase::P_PRECURSORS, on);
      }
      // features
      else if (action == dm_hulls_2d_)
      {
        win->canvas()->setLayerFlag(LayerDataBase::F_HULLS, on);
      }
      else if (action == dm_hull_2d_)
      {
        win->canvas()->setLayerFlag(LayerDataBase::F_HULL, on);
      }
      // consensus features
      else if (action == dm_elements_2d_)
      {
        win->canvas()->setLayerFlag(LayerDataBase::C_ELEMENTS, on);
      }
      // identifications
      else if (action == dm_ident_2d_)
      {
        win->canvas()->setLayerFlag(LayerDataBase::I_PEPTIDEMZ, on);
      }
    }
  }

  void TOPPViewBase::updateBarsAndMenus()
  {
    // Update filter bar, spectrum bar and layer bar
    layerActivated();
    updateMenu();
  }

  void TOPPViewBase::updateToolBar()
  {
    tool_bar_1d_->hide();
    tool_bar_2d_peak_->hide();
    tool_bar_2d_feat_->hide();
    tool_bar_2d_cons_->hide();
    tool_bar_2d_ident_->hide();

    PlotWidget* w = getActivePlotWidget();

    if (w)
    {
      // set intensity mode
      if (intensity_button_group_->button(w->canvas()->getIntensityMode()))
      {
        intensity_button_group_->button(w->canvas()->getIntensityMode())->setChecked(true);
      }
      else
      {
        log_->appendNewHeader(LogWindow::LogState::CRITICAL, OPENMS_PRETTY_FUNCTION, "Button for intensity mode does not exist");
      }
    }

    // 1D
    Plot1DWidget* w1 = getActive1DWidget();
    if (w1)
    {
      // draw mode
      draw_group_1d_->button(w1->canvas()->getDrawMode())->setChecked(true);

      // show/hide toolbars and buttons
      tool_bar_1d_->show();
    }

    // 2D
    Plot2DWidget* w2 = getActive2DWidget();
    if (w2)
    {
      // check if there is a layer before requesting data from it
      if (w2->canvas()->getLayerCount() > 0)
      {
        // peak draw modes
        if (w2->canvas()->getCurrentLayer().type == LayerDataBase::DT_PEAK)
        {
          dm_precursors_2d_->setChecked(w2->canvas()->getLayerFlag(LayerDataBase::P_PRECURSORS));
          tool_bar_2d_peak_->show();
        }
        // feature draw modes
        else if (w2->canvas()->getCurrentLayer().type == LayerDataBase::DT_FEATURE)
        {
          dm_hulls_2d_->setChecked(w2->canvas()->getLayerFlag(LayerDataBase::F_HULLS));
          dm_hull_2d_->setChecked(w2->canvas()->getLayerFlag(LayerDataBase::F_HULL));
          dm_unassigned_2d_->setChecked(w2->canvas()->getLayerFlag(LayerDataBase::F_UNASSIGNED));
          dm_label_2d_->setChecked(w2->canvas()->getCurrentLayer().label != LayerDataBase::L_NONE);
          tool_bar_2d_feat_->show();
        }
        // consensus feature draw modes
        else if (w2->canvas()->getCurrentLayer().type == LayerDataBase::DT_CONSENSUS)
        {
          dm_elements_2d_->setChecked(w2->canvas()->getLayerFlag(LayerDataBase::C_ELEMENTS));
          tool_bar_2d_cons_->show();
        }
        else if (w2->canvas()->getCurrentLayer().type == LayerDataBase::DT_IDENT)
        {
          dm_ident_2d_->setChecked(w2->canvas()->getLayerFlag(LayerDataBase::I_PEPTIDEMZ));
          tool_bar_2d_ident_->show();
        }
      }
    }

    // 3D
    Plot3DWidget* w3 = getActive3DWidget();
    if (w3)
    {
      // show no toolbars and buttons
    }
  }

  void TOPPViewBase::updateLayerBar()
  {
    layers_view_->update(getActivePlotWidget());
  }

  void TOPPViewBase::updateViewBar()
  {
    selection_view_->callUpdateEntries();
  }

  void TOPPViewBase::updateMenu()
  {
    FS_TV fs;
    LayerDataBase::DataType layer_type = LayerDataBase::DT_UNKNOWN;
    // is there a canvas?
    if (getActiveCanvas() != nullptr)
    {  
      fs |= TV_STATUS::HAS_CANVAS;
      // is there a layer?
      if (getActiveCanvas()->getLayerCount() != 0) 
      {
        fs |= TV_STATUS::HAS_LAYER;
        layer_type = getCurrentLayer()->type;
      }
    }
    // is this a 1D view
    if (getActive1DWidget() != nullptr) fs |= TV_STATUS::IS_1D_VIEW;
    // are we in 1D mirror mode
    if (getActive1DWidget() && getActive1DWidget()->canvas()->mirrorModeActive()) fs |= TV_STATUS::HAS_MIRROR_MODE;
    // is there a TOPP tool running
    if (topp_.process == nullptr)  fs |= TV_STATUS::TOPP_IDLE;
    
    menu_.update(fs, layer_type);
  }

  void TOPPViewBase::updateFilterBar()
  {
    PlotCanvas* canvas = getActiveCanvas();
    if (canvas == nullptr)
      return;

    if (canvas->getLayerCount() == 0)
      return;
    
    filter_list_->set(getActiveCanvas()->getCurrentLayer().filters);
  }

  void TOPPViewBase::layerFilterVisibilityChange(bool on) const
  {
    if (getActiveCanvas())
    {
      getActiveCanvas()->changeLayerFilterState(getActiveCanvas()->getCurrentLayerIndex(), on);
    }
  }

  void TOPPViewBase::layerActivated()
  {
    updateLayerBar();
    updateToolBar();
    updateViewBar();
    updateCurrentPath();
    updateFilterBar();
  }

  void TOPPViewBase::linkZoom()
  {
    zoom_together_ = !zoom_together_;
  }

  void TOPPViewBase::zoomOtherWindows() const
  {
    if (!zoom_together_) return;

    QList<QMdiSubWindow *> windows = ws_.subWindowList();
    if (!windows.count()) return;

    PlotWidget* w = getActivePlotWidget();
    auto new_visible_area = w->canvas()->getVisibleArea();
    // only zoom if other window is also (not) a chromatogram
    bool sender_is_chrom = w->canvas()->getCurrentLayer().type == LayerDataBase::DT_CHROMATOGRAM;

    // go through all windows, adjust the visible area where necessary
    for (int i = 0; i < int(windows.count()); ++i)
    {
      PlotWidget* specwidg = qobject_cast<PlotWidget*>(windows.at(i)->widget());
      if (!specwidg) continue;

      bool is_chrom = specwidg->canvas()->getCurrentLayer().type == LayerDataBase::DT_CHROMATOGRAM;
      if (is_chrom != sender_is_chrom) continue;
      // not the same dimensionality (e.g. Plot1DCanvas vs. 2DCanvas)
      if (w->canvas()->getName() != specwidg->canvas()->getName()) continue;

      specwidg->canvas()->setVisibleArea(new_visible_area);
    }
  }

  void TOPPViewBase::layerDeactivated()
  {
  }

  void TOPPViewBase::showPlotWidgetInWindow(PlotWidget* sw)
  {
    ws_.addSubWindow(sw);

    connect(sw->canvas(), &PlotCanvas::preferencesChange, this, &TOPPViewBase::updateLayerBar);
    connect(sw->canvas(), &PlotCanvas::layerActivated, this, &TOPPViewBase::layerActivated);
    connect(sw->canvas(), &PlotCanvas::layerModficationChange, this, &TOPPViewBase::updateLayerBar);
    connect(sw->canvas(), &PlotCanvas::layerZoomChanged, this, &TOPPViewBase::zoomOtherWindows);
    connect(sw, &PlotWidget::sendStatusMessage, this, &TOPPViewBase::showStatusMessage);
    connect(sw, &PlotWidget::sendCursorStatus, this, &TOPPViewBase::showCursorStatus);
    connect(sw, &PlotWidget::dropReceived, this, &TOPPViewBase::copyLayer);

    auto base_name = sw->canvas()->getCurrentLayer().getDecoratedName();

    // 1D spectrum specific signals
    Plot1DWidget* sw1 = qobject_cast<Plot1DWidget*>(sw);
    if (sw1)
    {
      connect(sw1, &Plot1DWidget::showCurrentPeaksAs2D, this, &TOPPViewBase::showCurrentPeaksAs2D);
      connect(sw1, &Plot1DWidget::showCurrentPeaksAs3D, this, &TOPPViewBase::showCurrentPeaksAs3D);
      connect(sw1, &Plot1DWidget::showCurrentPeaksAsIonMobility, this, &TOPPViewBase::showCurrentPeaksAsIonMobility);
      connect(sw1, &Plot1DWidget::showCurrentPeaksAsDIA, this, &TOPPViewBase::showCurrentPeaksAsDIA);
      base_name += " (1D)";
    }

    // 2D spectrum specific signals
    Plot2DWidget* sw2 = qobject_cast<Plot2DWidget*>(sw);
    if (sw2)
    {
      connect(sw2->getProjectionOntoX(), &Plot1DWidget::sendCursorStatus, this, &TOPPViewBase::showCursorStatus);
      connect(sw2->getProjectionOntoY(), &Plot1DWidget::sendCursorStatus, this, &TOPPViewBase::showCursorStatus);
      connect(sw2, &Plot2DWidget::showSpectrumAsNew1D, selection_view_, &DataSelectionTabs::showSpectrumAsNew1D);
      connect(sw2, &Plot2DWidget::showCurrentPeaksAsIonMobility, this, &TOPPViewBase::showCurrentPeaksAsIonMobility);
      connect(sw2, &Plot2DWidget::showCurrentPeaksAs3D, this, &TOPPViewBase::showCurrentPeaksAs3D);
      base_name += " (2D)";
    }

    // 3D spectrum specific signals
    Plot3DWidget* sw3 = qobject_cast<Plot3DWidget*>(sw);
    if (sw3)
    {
      connect(sw3, &Plot3DWidget::showCurrentPeaksAs2D,this, &TOPPViewBase::showCurrentPeaksAs2D);
      base_name += " (3D)";
    }

    sw->setWindowTitle(base_name.toQString());
    sw->addToTabBar(&tab_bar_, base_name, true);
    
    // show first window maximized (only visible windows are in the list)
    if (ws_.subWindowList().count() == 1)
    {
      sw->showMaximized();
    }
    else
    {
      sw->show();
    }
    showWindow(sw->getWindowId());
  }

  void TOPPViewBase::showGoToDialog() const
  {
    PlotWidget* w = getActivePlotWidget();
    if (w)
    {
      w->showGoToDialog();
    }
  }

  EnhancedWorkspace* TOPPViewBase::getWorkspace()
  {
    return &ws_;
  }

  PlotWidget* TOPPViewBase::getActivePlotWidget() const
  {
    // If the MDI window that holds all the subwindows for layers/spectra
    // is out-of-focus (e.g. because the table below was clicked and you moved out and into TOPPView),
    // currentSubWindow returns nullptr (i.e. no window is ACTIVE). In this case we get the one that is active
    // in the tabs (which SHOULD in theory be in-sync; due to a bug the way subwindow->tab does not work).
    // TODO check if we can reactivate automatically (e.g. double-check when TOPPView reacquires focus)
    if (!ws_.currentSubWindow())
    {
      // TODO think about using lastActivatedSubwindow_
      const int id = tab_bar_.currentIndex();

      if (id < 0 || id >= ws_.subWindowList().size()) return nullptr;

      return qobject_cast<PlotWidget*>(ws_.subWindowList()[id]->widget());
    }
    return qobject_cast<PlotWidget*>(ws_.currentSubWindow()->widget());
  }

  PlotCanvas* TOPPViewBase::getActiveCanvas() const
  {
    PlotWidget* sw = getActivePlotWidget();
    if (sw == nullptr)
    {
      return nullptr;
    }
    return sw->canvas();
  }

  Plot1DWidget* TOPPViewBase::getActive1DWidget() const
  {
    return qobject_cast<Plot1DWidget*>(getActivePlotWidget());
  }

  Plot2DWidget* TOPPViewBase::getActive2DWidget() const
  {
    return qobject_cast<Plot2DWidget*>(getActivePlotWidget());
  }

  Plot3DWidget* TOPPViewBase::getActive3DWidget() const
  {
    return qobject_cast<Plot3DWidget*>(getActivePlotWidget());
  }

  void TOPPViewBase::loadPreferences(String filename)
  {
    // compose default ini file path
    String default_ini_file = String(QDir::homePath()) + "/.TOPPView.ini";

    bool tool_params_added = false;

    if (filename == "") { filename = default_ini_file; }

    // load preferences, if file exists
    if (File::exists(filename))
    {
      bool error = false;
      Param tmp;
      try // the file might be corrupt
      {
        ParamXMLFile().load(filename, tmp);
      }
      catch (...)
      {
        error = true;
      }
      // apply preferences if they are of the current TOPPView version
      if (!error && tmp.exists("preferences:version") &&
          tmp.getValue("preferences:version").toString() == VersionInfo::getVersion())
      {
        try
        {
          setParameters(tmp.copy("preferences:"));
        }
        catch (Exception::InvalidParameter& /*e*/)
        {
          error = true;
        }
      }
      else
      {
        error = true;
      }

      // set parameters to defaults when something is fishy with the parameters file
      if (error)
      {
        // reset parameters (they will be stored again when TOPPView quits)
        setParameters(Param());
        cerr << "The TOPPView preferences files '" << filename << "' was ignored. It is no longer compatible with this TOPPView version and will be replaced." << endl;
      }
      else
      {
        // Load tool/util params
        if (scan_mode_ != TOOL_SCAN::FORCE_SCAN && tmp.hasSection("tool_params:"))
        {
          param_.insert("tool_params:", tmp.copy("tool_params:", true));
          tool_params_added = true;
        }
        // If the saved plugin path does not exist
        if (!tool_scanner_.setPluginPath(param_.getValue(user_section + "plugins_path").toString()))
        {
          // reset it to the default
          param_.setValue(user_section + "plugins_path", File::getUserDirectory() + "OpenMS_Plugins");
        }
      }
    }
    else if (filename != default_ini_file)
    {
      cerr << "Unable to load INI File: '" << filename << "'" << endl;
    }
    // Scan for tools if scan_mode is set to FORCE_SCAN or if the tool/util params could not be added for whatever reason
    if (!tool_params_added && scan_mode_ != TOOL_SCAN::SKIP_SCAN)
    {
      tool_scanner_.loadToolParams();
    }

    param_.setValue("PreferencesFile", filename);

    // set the recent files
    recent_files_.setFromParam(param_.copy("preferences:RecentFiles"));
  }

  void TOPPViewBase::savePreferences()
  {
    // replace recent files
    param_.removeAll("preferences:RecentFiles");
    param_.insert("preferences:RecentFiles:", recent_files_.getAsParam());

    // set version
    param_.setValue("preferences:version", VersionInfo::getVersion());
    // Make sure TOPP tool/util params have been inserted
    if (!param_.hasSection("tool_params:") && scan_mode_ != TOOL_SCAN::SKIP_SCAN)
    {
      tool_scanner_.waitForToolParams();
      param_.insert("tool_params:", tool_scanner_.getToolParams());
    }
    // check if the plugin path exists
    if (!tool_scanner_.setPluginPath(param_.getValue(user_section + "plugins_path").toString()))
    {
      // reset if it does not
      param_.setValue(user_section + "plugins_path", tool_scanner_.getPluginPath());
    }

    // save only the subsection that begins with "preferences:" and all tool params ("tool_params:")
    try
    {
      Param p;
      p.insert("preferences:", param_.copy("preferences:", true));
      p.insert("tool_params:", param_.copy("tool_params:", true));
      ParamXMLFile().store(string(param_.getValue("PreferencesFile")), p);
    }
    catch (Exception::UnableToCreateFile& /*e*/)
    {
      cerr << "Unable to create INI File: '" << string(param_.getValue("PreferencesFile")) << "'" << endl;
    }
  }

  QStringList TOPPViewBase::chooseFilesDialog_(const String& path_overwrite)
  {
    QString open_path = current_path_.toQString();
    if (!path_overwrite.empty())
    {
      open_path = path_overwrite.toQString();
    }
    // we use the QT file dialog instead of using QFileDialog::Names(...)
    // On Windows and Mac OS X, this static function will use the native file dialog and not a QFileDialog,
    // which prevents us from doing GUI testing on it.
    QFileDialog dialog(this, "Open file(s)", open_path, supported_types.toFileDialogFilter(FilterLayout::BOTH, true).toQString());
    dialog.setFileMode(QFileDialog::ExistingFiles);
    if (dialog.exec())
    {
       return dialog.selectedFiles();
    }
    return {};
  }

  void TOPPViewBase::openFilesByDialog(const String& dir)
  {
    for (const QString& filename : chooseFilesDialog_(dir))
    {
      addDataFile(filename, true, true);
    }
  }
  
  void TOPPViewBase::showTOPPDialog()
  {
    QAction* action = qobject_cast<QAction*>(sender());
    showTOPPDialog_(action->data().toBool());
  }

  void TOPPViewBase::showTOPPDialog_(bool visible_area_only)
  {
    // warn if hidden layer => wrong layer selected...
    const LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();
    if (!layer.visible)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "The current layer is not visible", "Have you selected the right layer for this action?");
    }

    // create and store unique file name prefix for files
    topp_.file_name = File::getTempDirectory() + "/TOPPView_" + File::getUniqueName();
    // Figure out the correct extension to give the temp file TODO start using OMS and cachedmzml

    if (!File::writable(topp_.file_name + "_ini"))
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Cannot create temporary file", String("Cannot write to '") + topp_.file_name + "'_ini!");
      return;
    }
    if (!param_.hasSection("tool_params:"))
    {
      tool_scanner_.waitForToolParams();
      param_.insert("tool_params:", tool_scanner_.getToolParams());
    }

    ToolsDialog tools_dialog(this, param_,
                             topp_.file_name + "_ini", current_path_, layer.type,
                             layer.getName(), &tool_scanner_);

    if (tools_dialog.exec() == QDialog::Accepted)
    {
      // Store tool name, input parameter and output parameter
      topp_.tool = tools_dialog.getTool();
      topp_.in = tools_dialog.getInput();
      topp_.out = tools_dialog.getOutput();
      topp_.visible_area_only = visible_area_only;
      // Build the input file name
      String file_extension;
      switch (layer.type)
        {
          case LayerDataBase::DataType::DT_PEAK:
            file_extension = FileTypes::typeToName(FileTypes::MZML);
            break;
          case LayerDataBase::DataType::DT_CHROMATOGRAM:
            file_extension = FileTypes::typeToName(FileTypes::MZML);
            break;
          case LayerDataBase::DataType::DT_FEATURE:
            file_extension = FileTypes::typeToName(FileTypes::FEATUREXML);
            break;
          case LayerDataBase::DataType::DT_CONSENSUS:
            file_extension = FileTypes::typeToName(FileTypes::CONSENSUSXML);
            break;
          case LayerDataBase::DataType::DT_IDENT:
            file_extension = FileTypes::typeToName(FileTypes::IDXML);
            break;
          default:
            file_extension = FileTypes::typeToName(FileTypes::UNKNOWN);
        }
      topp_.file_name_in = topp_.file_name + "_in." + file_extension;
      // Get the output file extension
      topp_.file_name_out = topp_.file_name + "_out." + tools_dialog.getExtension();
      // run the tool
      runTOPPTool_();
    }
  }

  void TOPPViewBase::rerunTOPPTool()
  {
    if (topp_.tool.empty())
    {
      QMessageBox::warning(this, "Error", "No TOPP tool was run before. Please run a tool first.");
      return;
    }
    // warn if hidden layer => wrong layer selected...
    const LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();
    if (!layer.visible)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "The current layer is not visible", "Have you selected the right layer for this action?");
    }

    // run the tool
    runTOPPTool_();
  }

  void TOPPViewBase::runTOPPTool_()
  {
    const LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();


    // delete old input and output file
    File::remove(topp_.file_name_in);
    File::remove(topp_.file_name_out);

    // test if files are writable
    if (!File::writable(topp_.file_name_in))
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Cannot create temporary file", String("Cannot write to '") + topp_.file_name_in + "'!");
      return;
    }
    if (!File::writable(topp_.file_name_out))
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Cannot create temporary file", String("Cannot write to '") + topp_.file_name_out + "'!");
      return;
    }

    // store data
    topp_.layer_name = layer.getName();
    topp_.window_id = getActivePlotWidget()->getWindowId();
    if (auto layer_1d = dynamic_cast<const LayerData1DBase*>(&layer))
    {
      topp_.spectrum_id = layer_1d->getCurrentIndex();
    }

    { // just a local scope
      auto visitor_data = topp_.visible_area_only
                          ? layer.storeVisibleData(getActiveCanvas()->getVisibleArea().getAreaUnit(), layer.filters)
                          : layer.storeFullData();
      visitor_data->saveToFile(topp_.file_name_in, ProgressLogger::GUI);
    }

    // compose argument list
    QStringList args;
    args << "-ini"
         << (topp_.file_name + "_ini").toQString()
         << QString("-%1").arg(topp_.in.toQString())
         << topp_.file_name_in.toQString()
         << "-no_progress";
    if (topp_.out != "")
    {
      args << QString("-%1").arg(topp_.out.toQString())
           << topp_.file_name_out.toQString();
    }

    // start log and show it
    log_->appendNewHeader(LogWindow::LogState::NOTICE, QString("Starting '%1'").arg(topp_.tool.toQString()), "");

    // initialize process
    topp_.process = new QProcess();
    topp_.process->setProcessChannelMode(QProcess::MergedChannels);

    // connect slots
    connect(topp_.process, &QProcess::readyReadStandardOutput, this, &TOPPViewBase::updateProcessLog);
    connect(topp_.process, CONNECTCAST(QProcess, finished, (int, QProcess::ExitStatus)), this, &TOPPViewBase::finishTOPPToolExecution);
    QString tool_executable = String(tool_scanner_.findPluginExecutable(topp_.tool)).toQString();
    if (tool_executable.isEmpty())
    {
      try
      {
        // find correct location of TOPP tool
        tool_executable = File::findSiblingTOPPExecutable(topp_.tool).toQString();
      }
      catch (Exception::FileNotFound & /*ex*/)
      {
        log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Could not locate executable!",
                              QString("Finding executable of TOPP tool '%1' failed. Please check your TOPP/OpenMS installation. Workaround: Add the bin/ directory to your PATH").arg(
                                      topp_.tool.toQString()));
        return;
      }
    }

    // update menu entries according to new state
    updateMenu();

    // start the actual process
    topp_.timer.restart();
    topp_.process->start(tool_executable, args);
    topp_.process->waitForStarted();

    if (topp_.process->error() == QProcess::FailedToStart)
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, QString("Failed to execute '%1'").arg(topp_.tool.toQString()), QString("Execution of TOPP tool '%1' failed with error: %2").arg(topp_.tool.toQString(), topp_.process->errorString()));

      // ensure that all tool output is emitted into log screen
      updateProcessLog();

      // re-enable Apply TOPP tool menus
      delete topp_.process;
      topp_.process = nullptr;
      updateMenu();
    }
  }

  void TOPPViewBase::finishTOPPToolExecution(int, QProcess::ExitStatus)
  {
    log_->addNewline();
    if (topp_.process->exitStatus() == QProcess::CrashExit)
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, QString("Execution of '%1' not successful!").arg(topp_.tool.toQString()),
                      QString("The tool crashed during execution. If you want to debug this crash, check the input files in '%1'"
                              " or enable 'debug' mode in the TOPP ini file.").arg(File::getTempDirectory().toQString()));
    }
    else if (topp_.process->exitCode() != 0) // NormalExit with non-zero exit code
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, QString("Execution of '%1' not successful!").arg(topp_.tool.toQString()),
                            QString("The tool ended with a non-zero exit code of '%1'. ").arg(topp_.process->exitCode()) +
                            QString("If you want to debug this, check the input files in '%1' or"
                                    " enable 'debug' mode in the TOPP ini file.").arg(File::getTempDirectory().toQString()));
    }
    else if (!topp_.out.empty())
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, QString("'%1' finished successfully").arg(topp_.tool.toQString()),
                      QString("Execution time: %1 ms").arg(topp_.timer.elapsed()));
      if (!File::readable(topp_.file_name_out))
      {
        log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Cannot read TOPP output", String("Cannot read '") + topp_.file_name_out + "'!");
      }
      else
      {
        addDataFile(topp_.file_name_out, true, false, topp_.layer_name + " (" + topp_.tool + ")", topp_.window_id, topp_.spectrum_id);
      }
    }

    // clean up
    delete topp_.process;
    topp_.process = nullptr;
    updateMenu();

    // clean up temporary files
    if (param_.getValue("preferences:topp_cleanup") == "true")
    {
      File::remove(topp_.file_name + "_ini");
      File::remove(topp_.file_name_in);
      File::remove(topp_.file_name_out);
    }
  }

  const LayerDataBase* TOPPViewBase::getCurrentLayer() const
  {
    PlotCanvas* canvas = getActiveCanvas();
    if (canvas == nullptr)
    {
      return nullptr;
    }
    return &(canvas->getCurrentLayer());
  }

  LayerDataBase* TOPPViewBase::getCurrentLayer()
  {
    PlotCanvas* canvas = getActiveCanvas();
    if (canvas == nullptr)
    {
      return nullptr;
    }
    return &(canvas->getCurrentLayer());
  }

  void TOPPViewBase::toggleProjections()
  {
    Plot2DWidget* w = getActive2DWidget();
    if (w)
    {
      // update minimum size before
      if (!w->projectionsVisible())
      {
        setMinimumSize(700, 700);
      }
      else
      {
        setMinimumSize(400, 400);
      }
      w->toggleProjections();
    }
  }

  void TOPPViewBase::annotateWithAMS()
  { // this should only be callable if current layer's type is of type DT_PEAK
    LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();
    LayerAnnotatorAMS annotator(this);
    assert(log_ != nullptr);
    if (!annotator.annotateWithFileDialog(layer, *log_, current_path_))
    {
      return;
    }
  }

  void TOPPViewBase::annotateWithID()
  { // this should only be callable if current layer's type is one of DT_PEAK, DT_FEATURE, DT_CONSENSUS
    LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();
    LayerAnnotatorPeptideID annotator(this);
    assert(log_ != nullptr);
    if (!annotator.annotateWithFileDialog(layer, *log_, current_path_))
    {
      return;
    }
    selection_view_->setCurrentIndex(DataSelectionTabs::IDENT_IDX); //switch to ID view
    selection_view_->currentTabChanged(DataSelectionTabs::IDENT_IDX);
  }

  void TOPPViewBase::annotateWithOSW()
  { // this should only be callable if current layer's type is of type DT_CHROMATOGRAM
    LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();
    LayerAnnotatorOSW annotator(this);
    assert(log_ != nullptr);
    if (!annotator.annotateWithFileDialog(layer, *log_, current_path_))
    {
      return;
    }
    selection_view_->setCurrentIndex(DataSelectionTabs::DIAOSW_IDX); // switch to DIA view
    selection_view_->currentTabChanged(DataSelectionTabs::DIAOSW_IDX);
  }

  void TOPPViewBase::showSpectrumGenerationDialog()
  {
    // TheoreticalSpectrumGenerationDialog spec_gen_dialog;
    if (spec_gen_dialog_.exec())
    {
      // spectrum is generated in the dialog, so just receive it here
      PeakSpectrum spectrum = spec_gen_dialog_.getSpectrum();

      PeakMap new_exp;
      new_exp.addSpectrum(spectrum);
      ExperimentSharedPtrType new_exp_sptr(new PeakMap(new_exp));
      FeatureMapSharedPtrType f_dummy(new FeatureMapType());
      ConsensusMapSharedPtrType c_dummy(new ConsensusMapType());
      ODExperimentSharedPtrType od_dummy(new OnDiscMSExperiment());
      vector<PeptideIdentification> p_dummy;
      addData(f_dummy, c_dummy, p_dummy, new_exp_sptr, od_dummy, LayerDataBase::DT_PEAK, false, true, true, "", spec_gen_dialog_.getSequence() + " (theoretical)");

      // ensure spectrum is drawn as sticks
      draw_group_1d_->button(Plot1DCanvas::DM_PEAKS)->setChecked(true);
      setDrawMode1D(Plot1DCanvas::DM_PEAKS);
    }
  }

  void TOPPViewBase::showSpectrumAlignmentDialog()
  {
    Plot1DWidget* active_1d_window = getActive1DWidget();
    if (!active_1d_window || !active_1d_window->canvas()->mirrorModeActive())
    {
      return;
    }
    Plot1DCanvas* cc = active_1d_window->canvas();

    SpectrumAlignmentDialog spec_align_dialog(active_1d_window);
    if (spec_align_dialog.exec())
    {
      Int layer_index_1 = spec_align_dialog.get1stLayerIndex();
      Int layer_index_2 = spec_align_dialog.get2ndLayerIndex();

      // two layers must be selected:
      if (layer_index_1 < 0 || layer_index_2 < 0)
      {
        QMessageBox::information(this, "Layer selection invalid", "You must select two layers for an alignment.");
        return;
      }

      Param param;
      double tolerance = spec_align_dialog.getTolerance();
      param.setValue("tolerance", tolerance, "Defines the absolute (in Da) or relative (in ppm) mass tolerance");
      String unit_is_ppm = spec_align_dialog.isPPM() ? "true" : "false";
      param.setValue("is_relative_tolerance", unit_is_ppm, "If true, the mass tolerance is interpreted as ppm value otherwise in Dalton");

      active_1d_window->performAlignment((UInt)layer_index_1, (UInt)layer_index_2, param);

      double al_score = cc->getAlignmentScore();
      Size al_size = cc->getAlignmentSize();

      QMessageBox::information(this, "Alignment performed", QString("Aligned %1 pairs of peaks (Score: %2).").arg(al_size).arg(al_score));
    }
  }
  
  void TOPPViewBase::showCurrentPeaksAs2D()
  {
    LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();
    auto* lp = dynamic_cast<LayerDataPeak*>(&layer);
    if (!lp) return;

    ExperimentSharedPtrType exp_sptr = lp->getPeakDataMuteable();
    ODExperimentSharedPtrType od_exp_sptr = lp->getOnDiscPeakData();

    // open new 2D widget
    Plot2DWidget* w = new Plot2DWidget(getCanvasParameters(2), &ws_);

    // add data
    if (!w->canvas()->addPeakLayer(exp_sptr, od_exp_sptr, layer.filename))
    {
      return;
    }

    showPlotWidgetInWindow(w);
    updateMenu();
  }


  void TOPPViewBase::showCurrentPeaksAsIonMobility(const MSSpectrum& spec)
  {
    const LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();
    
    ExperimentSharedPtrType exp(new MSExperiment(IMDataConverter::reshapeIMFrameToMany(spec)));
    // hack, but currently not avoidable, because 2D widget does not support IM natively yet...
    // for (auto& spec : exp->getSpectra()) spec.setRT(spec.getDriftTime());

    // open new 2D widget
    Plot2DWidget* w = new Plot2DWidget(getCanvasParameters(2), &ws_);
    // map to IM + MZ
    w->setMapper(DimMapper<2>({IMTypes::fromIMUnit(exp->getSpectra()[0].getDriftTimeUnit()), DIM_UNIT::MZ}));

    // add data
    if (!w->canvas()->addPeakLayer(exp, PlotCanvas::ODExperimentSharedPtrType(new OnDiscMSExperiment()), layer.filename + " (IM Frame)"))
    {
      return;
    }

    showPlotWidgetInWindow(w);
    updateMenu();
  }

  void TOPPViewBase::showCurrentPeaksAsDIA(const Precursor& pc, const MSExperiment& exp)
  {
    const LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();
    auto* lp = dynamic_cast<const LayerDataPeak*>(&layer);
    if (!lp || !lp->isDIAData())
    {
      std::cout << "Layer does not contain DIA / SWATH-MS data" << std::endl;
      return;
    }

    // Add spectra into a MSExperiment, sort and prepare it for display
    ExperimentSharedPtrType tmpe(new OpenMS::MSExperiment() );

    // Collect all MS2 spectra with the same precursor as the current spectrum
    // (they are in the same SWATH window)
    String caption_add = "";

    double lower = pc.getMZ() - pc.getIsolationWindowLowerOffset();
    double upper = pc.getMZ() + pc.getIsolationWindowUpperOffset();

    Size k = 0;
    for (const auto& spec : exp)
    {
      if (spec.getMSLevel() == 2 && !spec.getPrecursors().empty() )
      {
        if (fabs(spec.getPrecursors()[0].getMZ() - pc.getMZ()) < 1e-4)
        {
          // Get the spectrum in question (from memory or disk) and add to
          // the newly created MSExperiment
          if (spec.size() > 0)
          {
            // Get data from memory - copy data and tell TOPPView that this
            // is MS1 data so that it will be displayed properly in 2D and 3D
            // view
            MSSpectrum t = spec;
            t.setMSLevel(1);
            tmpe->addSpectrum(t);
          }
          else if (lp->getOnDiscPeakData()->getNrSpectra() > k)
          {
            // Get data from disk - copy data and tell TOPPView that this is
            // MS1 data so that it will be displayed properly in 2D and 3D
            // view
            MSSpectrum t = lp->getOnDiscPeakData()->getSpectrum(k);
            t.setMSLevel(1);
            tmpe->addSpectrum(t);
          }
        }
      }
      k++;
    }
    caption_add = "(DIA window " + String(lower) + " - " + String(upper) + ")";
    
    tmpe->sortSpectra();
    tmpe->updateRanges();

    // open new 2D widget
    Plot2DWidget* w = new Plot2DWidget(getCanvasParameters(2), &ws_);

    // add data
    if (!w->canvas()->addPeakLayer(tmpe, PlotCanvas::ODExperimentSharedPtrType(new OnDiscMSExperiment()), layer.filename))
    {
      return;
    }

    String caption = layer.getName();
    caption += caption_add;
    w->canvas()->setLayerName(w->canvas()->getCurrentLayerIndex(), caption);
    showPlotWidgetInWindow(w);
    updateMenu();
  }

  void TOPPViewBase::showCurrentPeaksAs3D()
  {
    // we first pick the layer with 3D support which is closest (or ideally identical) to the currently active layer
    // we might find that there is no compatible layer though...
    // if some day more than one type of data is supported, we need to adapt the code below accordingly

    const int BIGINDEX = 10000; // some large number which is never reached
    const int target_layer = (int) getActiveCanvas()->getCurrentLayerIndex();
    int best_candidate = BIGINDEX;
    for (int i = 0; i < (int) getActiveCanvas()->getLayerCount(); ++i)
    {
      if ((LayerDataBase::DT_PEAK == getActiveCanvas()->getLayer(i).type) && // supported type
          (std::abs(i - target_layer) < std::abs(best_candidate - target_layer))) // & smallest distance to active layer
      {
        best_candidate = i;
      }
    }

    if (best_candidate == BIGINDEX)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "No compatible layer",
          "No layer found which is supported by the 3D view.");
      return;
    }


    if (best_candidate != target_layer)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "Auto-selected compatible layer",
          "The currently active layer cannot be viewed in 3D view. The closest layer which is supported by the 3D view was selected!");
    }

    LayerDataBase& layer = const_cast<LayerDataBase&>(getActiveCanvas()->getLayer(best_candidate));
    auto* lp = dynamic_cast<LayerDataPeak*>(&layer);
    if (!lp)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "Wrong layer type", "Something went wrong during layer selection. Please report this problem with a description of your current layers!");
      return;
    }
    // open new 3D widget
    Plot3DWidget* w = new Plot3DWidget(getCanvasParameters(3), &ws_);

    ExperimentSharedPtrType exp_sptr = lp->getPeakDataMuteable();

    if (lp->isIonMobilityData())
    {
      // Determine ion mobility unit (default is milliseconds)
      String unit = "ms";
      if (exp_sptr->metaValueExists("ion_mobility_unit"))
      {
        unit = exp_sptr->getMetaValue("ion_mobility_unit");
      }
      String label = "Ion Mobility [" + unit + "]";

      w->canvas()->openglwidget()->setYLabel(label.c_str());
    }

    if (!w->canvas()->addPeakLayer(exp_sptr, PlotCanvas::ODExperimentSharedPtrType(new OnDiscMSExperiment()), layer.filename))
    {
      return;
    }

    if (getActive1DWidget()) // switch from 1D to 3D
    {
      // TODO:
      //- doesn't make sense for fragment scan
      //- build new Area with mz range equal to 1D visible range
      //- rt range either overall MS1 data range or some convenient window

    }
    else if (getActive2DWidget()) // switch from 2D to 3D
    {
      w->canvas()->setVisibleArea(getActiveCanvas()->getVisibleArea());
    }

    showPlotWidgetInWindow(w);

    // set intensity mode (after spectrum has been added!)
    setIntensityMode(PlotCanvas::IM_SNAP);
    updateMenu();
  }

  void TOPPViewBase::updateProcessLog()
  {
    log_->appendText(topp_.process->readAllStandardOutput());
  }

  Param TOPPViewBase::getCanvasParameters(UInt dim) const
  {
    Param out = param_.copy(String(user_section + "") + dim + "d:", true);
    out.setValue("default_path", param_.getValue(user_section + "default_path").toString());
    return out;
  }

  void TOPPViewBase::abortTOPPTool()
  {
    if (topp_.process)
    {
      // block signals to avoid error message from finished() signal
      topp_.process->blockSignals(true);
      // kill and delete the process
      topp_.process->terminate();
      delete topp_.process;
      topp_.process = nullptr;

      // finish log with new line
      log_->addNewline();

      updateMenu();
    }
  }
  
  void TOPPViewBase::loadFiles(const StringList& list, QSplashScreen* splash_screen)
  {
    static StringList colors = { "@bw", "@bg", "@b", "@r", "@g", "@m" };
    static StringList gradients = { "Linear|0,#ffffff;100,#000000" , "Linear|0,#dddddd;100,#000000" , "Linear|0,#000000;100,#000000",
                                    "Linear|0,#ff0000;100,#ff0000" , "Linear|0,#00ff00;100,#00ff00" , "Linear|0,#ff00ff;100,#ff00ff" };
    bool last_was_plus = false;
    bool last_was_annotation = false;
    for (StringList::const_iterator it = list.begin(); it != list.end(); ++it)
    {
      if (*it == "+")
      {
        last_was_plus = true;
        continue;
      }
      if (*it == "!")
      {
        last_was_annotation = true;
        continue;
      }

      // no matter what the current item is, after we are done with it, 
      // we need to reset the 'glue' symbols
      RAIICleanup reset([&]() {
        last_was_plus = false;
        last_was_annotation = false;
      });

      if (std::find(colors.begin(), colors.end(), *it) != colors.end())
      { // its a color!
        if ((getActive2DWidget() != nullptr || getActive3DWidget() != nullptr) && getActiveCanvas() != nullptr)
        {
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
          tmp.setValue("dot:gradient", gradients[Helpers::indexOf(colors, *it)]);
          getActiveCanvas()->setCurrentLayerParameters(tmp);
        }
        continue;
      }
      
      splash_screen->showMessage((String("Loading file: ") + *it).toQString());
      splash_screen->repaint();
      QApplication::processEvents();

      if (!getActivePlotWidget())
      {
        if (last_was_annotation)
        {
          log_->appendNewHeader(LogWindow::LogState::WARNING, "Error", "Cannot annotate without having added layers before.");
          continue;
        }
        // create new tab (also in case of last_was_plus)...
        addDataFile(*it, false, true); // add data file but don't show options
        continue;
      }

      // we have an active widget
      if (last_was_plus)
      { // add to current tab
        addDataFile(*it, false, true, "", getActivePlotWidget()->getWindowId());
        continue;
      }
      else if (last_was_annotation)
      { // try to treat file as annotation file and annotate current layer
        auto l = getCurrentLayer();
        if (l)
        {
          auto annotator = LayerAnnotatorBase::getAnnotatorWhichSupports(*it);
          if (annotator.get() == nullptr)
          {
            log_->appendNewHeader(LogWindow::LogState::NOTICE, "Error", String("Filename '" + *it + "' has unsupported file type. No annotation performed.").toQString());
          }
          else
          { // we have an annotator ...
            annotator->annotateWithFilename(*l, *log_, *it); // ID tabs are automatically enabled
          }
        }
      }        
      else
      { // create new tab
        addDataFile(*it, false, true); // add data file but don't show options
      }
    }
  }

  void TOPPViewBase::saveLayerAll() const
  {
    getActiveCanvas()->saveCurrentLayer(false);
  }

  void TOPPViewBase::saveLayerVisible() const
  {
    getActiveCanvas()->saveCurrentLayer(true);
  }

  void TOPPViewBase::toggleGridLines() const
  {
    getActiveCanvas()->showGridLines(!getActiveCanvas()->gridLinesShown());
  }

  void TOPPViewBase::toggleAxisLegends() const
  {
    getActivePlotWidget()->showLegend(!getActivePlotWidget()->isLegendShown());
  }

  void TOPPViewBase::toggleInterestingMZs() const
  {
    auto w = getActive1DWidget();
    if (w == nullptr) return;
    w->canvas()->setDrawInterestingMZs(!w->canvas()->isDrawInterestingMZs());
  }

  void TOPPViewBase::showPreferences() const
  {
    getActiveCanvas()->showCurrentLayerPreferences();
  }

  void TOPPViewBase::metadataFileDialog()
  {
    QStringList files = chooseFilesDialog_();
    FileHandler fh;
    fh.getOptions().setMetadataOnly(true);
    for (QStringList::iterator it = files.begin(); it != files.end(); ++it)
    {
      ExperimentType exp;
      try
      {
        QMessageBox::critical(this, "Error", "Only raw data files (mzML, DTA etc) are supported to view their meta data.");
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::critical(this, "Error", (String("Error while reading data: ") + e.what()).c_str());
        return;
      }
      MetaDataBrowser dlg(false, this);
      dlg.add(exp);
      dlg.exec();
    }
  }


  void TOPPViewBase::showSpectrumMetaData(int spectrum_index) const
  {
    getActiveCanvas()->showMetaData(true, spectrum_index);
  }

  void TOPPViewBase::copyLayer(const QMimeData* data, QWidget* source, int id)
  {
    SpectraTreeTab* spec_view = (source ? qobject_cast<SpectraTreeTab*>(source->parentWidget()) : nullptr);
    try
    {
      // NOT USED RIGHT NOW, BUT KEEP THIS CODE (it was hard to find out how this is done)
      // decode data to get the row
      // QByteArray encoded_data = data->data(data->formats()[0]);
      // QDataStream stream(&encoded_data, QIODevice::ReadOnly);
      // int row, col;
      // stream >> row >> col;

      // set wait cursor
      setCursor(Qt::WaitCursor);
      RAIICleanup cl([&]() { setCursor(Qt::ArrowCursor); });

      // determine where to copy the data
      UInt new_id = (id == -1) ? 0 : id;

      if (source == layers_view_)
      {
        // only the selected row can be dragged => the source layer is the selected layer
        LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();

        // attach feature, consensus and peak data          (new OnDiscMSExperiment());
        FeatureMapSharedPtrType features(new FeatureMapType());
        if (auto* lp = dynamic_cast<LayerDataFeature*>(&layer)) features = lp->getFeatureMap();

        ConsensusMapSharedPtrType consensus(new ConsensusMapType());
        if (auto* lp = dynamic_cast<LayerDataConsensus*>(&layer)) consensus = lp->getConsensusMap();

        ExperimentSharedPtrType peaks(new ExperimentType());
        ODExperimentSharedPtrType on_disc_peaks(new OnDiscMSExperiment());
        if (auto* lp = dynamic_cast<LayerDataPeak*>(&layer))
        {
          peaks = lp->getPeakDataMuteable();
          on_disc_peaks = lp->getOnDiscPeakData();
        }
        if (auto* lp = dynamic_cast<LayerDataChrom*>(&layer))
        {
          peaks = lp->getChromatogramData();
          on_disc_peaks = lp->getOnDiscPeakData();
        }
        // if the layer provides identification data -> retrieve it
        vector<PeptideIdentification> peptides;
        if (auto p = dynamic_cast<IPeptideIds*>(&layer); p != nullptr)
        {
          peptides = p->getPeptideIds();
        }

        // add the data
        addData(features, consensus, peptides, peaks, on_disc_peaks, layer.type, false, false, true, layer.filename, layer.getName(), new_id);
      }
      else if (spec_view != nullptr)
      {
        ExperimentSharedPtrType new_exp_sptr(new ExperimentType());
        if (LayerDataBase::DataType current_type; spec_view->getSelectedScan(*new_exp_sptr, current_type))
        {
          ODExperimentSharedPtrType od_dummy(new OnDiscMSExperiment());
          FeatureMapSharedPtrType f_dummy(new FeatureMapType());
          ConsensusMapSharedPtrType c_dummy(new ConsensusMapType());
          vector<PeptideIdentification> p_dummy;
          const LayerDataBase& layer = getActiveCanvas()->getCurrentLayer();
          addData(f_dummy, c_dummy, p_dummy, new_exp_sptr, od_dummy, current_type, false, false, true, layer.filename, layer.getName(), new_id);
        }
      }
      else if (source == nullptr)
      {
        // drag source is external
        if (data->hasUrls())
        {
          QList<QUrl> urls = data->urls();
          // use a QTimer for external sources to make the source (e.g. Windows Explorer responsive again)
          // Using a QueuedConnection for the DragEvent does not solve the problem (Qt 5.15) -- see previous (reverted) commit
          QTimer::singleShot(50, [this, urls, new_id]() {
            for (const QUrl& url : urls)
            {
              addDataFile(url.toLocalFile(), false, true, "", new_id);
            }
          });
        }
      }
    }
    catch (Exception::BaseException& e)
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Error while creating layer", e.what());
    }
  }

  void TOPPViewBase::updateCurrentPath()
  {
    // do not update if the user disabled this feature.
    if (param_.getValue(user_section + "default_path_current") != "true")
    {
      return;
    }

    // reset
    current_path_ = param_.getValue(user_section + "default_path").toString();

    // update if the current layer has a path associated
    if (getActiveCanvas() && getActiveCanvas()->getLayerCount() != 0 && getActiveCanvas()->getCurrentLayer().filename != "")
    {
      current_path_ = File::path(getActiveCanvas()->getCurrentLayer().filename);
    }
  }

  void TOPPViewBase::showSpectrumBrowser()
  {
    views_dockwidget_->show();
    updateViewBar();
  }

  void TOPPViewBase::fileChanged_(const String& filename)
  {
    // check if file has been deleted
    if (!QFileInfo(filename.toQString()).exists())
    {
      watcher_->removeFile(filename);
      return;
    }

    // iterate over all windows and determine which need an update
    std::vector<std::pair<const PlotWidget*, Size> > needs_update;
    for (const auto& mdi_window : ws_.subWindowList())
    {
      const PlotWidget* sw = qobject_cast<const PlotWidget*>(mdi_window);
      if (sw == nullptr) return;

      Size lc = sw->canvas()->getLayerCount();
      // determine if widget stores one or more layers for the given filename (->needs update)
      for (Size j = 0; j != lc; ++j)
      {
        if (sw->canvas()->getLayer(j).filename == filename)
        {
          needs_update.emplace_back(sw, j);
        }
      }
    }

    if (needs_update.empty()) // no layer references data of filename
    {
      watcher_->removeFile(filename); // remove watcher
      return;
    }
    
    // std::cout << "Number of Layers that need update: " << needs_update.size() << std::endl;
    pair<const PlotWidget*, Size>& slp = needs_update[0];
    const PlotWidget* sw = slp.first;
    Size layer_index = slp.second;

    bool user_wants_update = false;
    if (param_.getValue(user_section + "on_file_change") == "update automatically") //automatically update
    {
      user_wants_update = true;
    }
    else if (param_.getValue(user_section + "on_file_change") == "ask") // ask the user if the layer should be updated
    {
      if (watcher_msgbox_ == true) // we already have a dialog for that opened... do not ask again
      {
        return;
      }
      // track that we will show the msgbox and we do not need to show it again if file changes once more and the dialog is still open
      watcher_msgbox_ = true;
      QMessageBox msg_box;
      QAbstractButton* ok = msg_box.addButton(QMessageBox::Ok);
      msg_box.addButton(QMessageBox::Cancel);
      msg_box.setWindowTitle("Layer data changed");
      msg_box.setText((String("The data of file '") + filename + "' has changed.<BR>Update layers?").toQString());
      msg_box.exec();
      watcher_msgbox_ = false;
      if (msg_box.clickedButton() == ok)
      {
        user_wants_update = true;
      }
    }

    if (user_wants_update == false)
    {
      return;
    }
    LayerDataBase& layer = sw->canvas()->getLayer(layer_index);
    // reload data
    if (auto* lp = dynamic_cast<LayerDataPeak*>(&layer)) // peak data
    {
      try
      {
        FileHandler().loadExperiment(layer.filename, *lp->getPeakDataMuteable(), {}, ProgressLogger::NONE, true, true);
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::critical(this, "Error", (String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
        lp->getPeakDataMuteable()->clear(true);
      }
      lp->getPeakDataMuteable()->sortSpectra(true);
      lp->getPeakDataMuteable()->updateRanges(1);
    }
    else if (auto* lp = dynamic_cast<LayerDataFeature*>(&layer)) // feature data
    {
      try
      {
        FileHandler().loadFeatures(layer.filename, *lp->getFeatureMap());
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::critical(this, "Error", (String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
        lp->getFeatureMap()->clear(true);
      }
      lp->getFeatureMap()->updateRanges();
    }
    else if (auto* lp = dynamic_cast<LayerDataConsensus*>(&layer)) // consensus feature data
    {
      try
      {
        FileHandler().loadConsensusFeatures(layer.filename, *lp->getConsensusMap(), {FileTypes::CONSENSUSXML});
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::critical(this, "Error", (String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
        lp->getConsensusMap()->clear(true);
      }
      lp->getConsensusMap()->updateRanges();
    }
    else if (auto* lp = dynamic_cast<LayerDataChrom*>(&layer)) // chromatogram
    {
      // TODO CHROM
      try
      {
        FileHandler().loadExperiment(layer.filename, *lp->getChromatogramData(), {}, ProgressLogger::NONE, true, true);
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::critical(this, "Error", (String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
        lp->getChromatogramData()->clear(true);
      }
      lp->getChromatogramData()->sortChromatograms(true);
      lp->getChromatogramData()->updateRanges(1);
    }

    // update all layers that need an update
    for (Size i = 0; i != needs_update.size(); ++i)
    {
      sw->canvas()->updateLayer(needs_update[i].second);
    }
    layerActivated();

    // temporarily remove and read filename from watcher_ as a workaround for bug #233
    // This might not be a 'bug' but rather unfortunate behaviour (even in Qt5) if the file was actually deleted and recreated by an external tool
    // (some TextEditors seem to do this), see https://stackoverflow.com/a/30076119;
    watcher_->removeFile(filename);
    watcher_->addFile(filename);
  }

  TOPPViewBase::~TOPPViewBase()
  {
    savePreferences();
    abortTOPPTool();
  }

} // namespace OpenMS
