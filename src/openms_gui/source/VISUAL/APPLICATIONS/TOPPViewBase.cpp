// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/EnumHelpers.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/DIALOGS/SpectrumAlignmentDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewPrefDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>
#include <OpenMS/VISUAL/LayerListView.h>
#include <OpenMS/VISUAL/LogWindow.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/VISUAL/SpectraIdentificationViewWidget.h>
#include <OpenMS/VISUAL/SpectraSelectionTabs.h>
#include <OpenMS/VISUAL/SpectraViewWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DWidget.h>

//Qt
#include <QtCore/QSettings>
#include <QtCore/QDate>
#include <QtCore/QDir>
#include <QtCore/QTime>
#include <QtCore/QUrl>
#include <QtWidgets/QCheckBox>
#include <QCloseEvent>
#include <QtWidgets/QDesktopWidget>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QInputDialog>
#include <QtWidgets/QMessageBox>
#include <QPainter>
#include <QtWidgets/QSplashScreen>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QToolTip>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QTreeWidget>
#include <QtWidgets/QTreeWidgetItem>
#include <QtWidgets/QWhatsThis>
#include <QTextCodec>

#include <boost/math/special_functions/fpclassify.hpp>

#include <algorithm>
#include <utility>

using namespace std;

namespace OpenMS
{
  using namespace Internal;
  using namespace Math;

  const String TOPPViewBase::CAPTION_3D_SUFFIX_ = " (3D)";

  FileTypes::FileTypeList<11> supported_types({ FileTypes::MZML, FileTypes::MZXML, FileTypes::MZDATA, FileTypes::SQMASS,
                                               FileTypes::FEATUREXML, FileTypes::CONSENSUSXML, FileTypes::IDXML,
                                               FileTypes::DTA, FileTypes::DTA2D,
                                               FileTypes::BZ2, FileTypes::GZ });
  TOPPViewBase::TOPPViewBase(QWidget* parent) :
    QMainWindow(parent),
    DefaultParamHandler("TOPPViewBase"),
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
    QRect screen_geometry = QApplication::desktop()->screenGeometry();
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

    connect(&ws_, &EnhancedWorkspace::subWindowActivated, [this](QMdiSubWindow* window) {
      if (window != nullptr) /* 0 upon terminate */ updateBarsAndMenus(); 
    });
    connect(&ws_, &EnhancedWorkspace::dropReceived, this, &TOPPViewBase::copyLayer);
    box_layout->addWidget(&ws_);

    //################## STATUS #################
    // create status bar
    message_label_ = new QLabel(statusBar());
    statusBar()->addWidget(message_label_, 1);

    rt_label_ = new QLabel("RT: 12345678", statusBar());
    rt_label_->setMinimumSize(rt_label_->sizeHint());
    rt_label_->setText("");
    statusBar()->addPermanentWidget(rt_label_, 0);
    mz_label_ = new QLabel("m/z: 123456780912", statusBar());
    mz_label_->setMinimumSize(mz_label_->sizeHint());
    mz_label_->setText("");
    statusBar()->addPermanentWidget(mz_label_, 0);

    //################## TOOLBARS #################
    //create toolbars and connect signals
    QToolButton* b;

    //--Basic tool bar for all views--
    tool_bar_ = addToolBar("Basic tool bar");
    tool_bar_->setObjectName("tool_bar");

    //intensity modes
    intensity_button_group_ = new QButtonGroup(tool_bar_);
    intensity_button_group_->setExclusive(true);

    b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/lin.png"));
    b->setToolTip("Intensity: Normal");
    b->setShortcut(Qt::Key_N);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Normal<BR><BR>Intensity is displayed unmodified.<BR>(Hotkey: N)");
    intensity_button_group_->addButton(b, SpectrumCanvas::IM_NONE);
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
    intensity_button_group_->addButton(b, SpectrumCanvas::IM_PERCENTAGE);
    tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/snap.png"));
    b->setToolTip("Intensity: Snap to maximum displayed intensity");
    b->setShortcut(Qt::Key_S);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Snap to maximum displayed intensity<BR><BR> In this mode the"
                    " color gradient is adapted to the maximum currently displayed intensity."
                    "<BR>(Hotkey: S)");
    intensity_button_group_->addButton(b, SpectrumCanvas::IM_SNAP);
    tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/log.png"));
    b->setToolTip("Intensity: Use log scaling for colors");
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Logarithmic scaling of intensities for color calculation");
    intensity_button_group_->addButton(b, SpectrumCanvas::IM_LOG);
    tool_bar_->addWidget(b);

    connect(intensity_button_group_, CONNECTCAST(QButtonGroup,buttonClicked,(int)), this, &TOPPViewBase::setIntensityMode);
    tool_bar_->addSeparator();

    //common buttons
    QAction* reset_zoom_button = tool_bar_->addAction(QIcon(":/reset_zoom.png"), "Reset Zoom", this, &TOPPViewBase::resetZoom);
    reset_zoom_button->setWhatsThis("Reset zoom: Zooms out as far as possible and resets the zoom history.<BR>(Hotkey: Backspace)");

    tool_bar_->show();

    //--1D toolbar--
    tool_bar_1d_ = addToolBar("1D tool bar");
    tool_bar_1d_->setObjectName("1d_tool_bar");

    //draw modes 1D
    draw_group_1d_ = new QButtonGroup(tool_bar_1d_);
    draw_group_1d_->setExclusive(true);

    b = new QToolButton(tool_bar_1d_);
    b->setIcon(QIcon(":/peaks.png"));
    b->setToolTip("Peak mode");
    b->setShortcut(Qt::Key_I);
    b->setCheckable(true);
    b->setWhatsThis("1D Draw mode: Peaks<BR><BR>Peaks are displayed as sticks.");
    draw_group_1d_->addButton(b, Spectrum1DCanvas::DM_PEAKS);
    tool_bar_1d_->addWidget(b);

    b = new QToolButton(tool_bar_1d_);
    b->setIcon(QIcon(":/lines.png"));
    b->setToolTip("Raw data mode");
    b->setShortcut(Qt::Key_R);
    b->setCheckable(true);
    b->setWhatsThis("1D Draw mode: Raw data<BR><BR>Peaks are displayed as a continuous line.");
    draw_group_1d_->addButton(b, Spectrum1DCanvas::DM_CONNECTEDLINES);
    tool_bar_1d_->addWidget(b);

    connect(draw_group_1d_, CONNECTCAST(QButtonGroup, buttonClicked, (int)), this, &TOPPViewBase::setDrawMode1D);
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
    //button menu
    group_label_2d_ = new QActionGroup(dm_label_2d_);
    QMenu* menu = new QMenu(dm_label_2d_);
    for (Size i = 0; i < LayerData::SIZE_OF_LABEL_TYPE; ++i)
    {
      QAction* temp = group_label_2d_->addAction(
        QString(LayerData::NamesOfLabelType[i].c_str()));
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
    //button menu
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
    selection_view_ = new SpectraSelectionTabs(views_dockwidget_, this);
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
    current_path_ = param_.getValue("preferences:default_path");

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
    //general
    defaults_.setValue("preferences:default_map_view", "2d", "Default visualization mode for maps.");
    defaults_.setValidStrings("preferences:default_map_view", ListUtils::create<String>("2d,3d"));
    defaults_.setValue("preferences:default_path", ".", "Default path for loading and storing files.");
    defaults_.setValue("preferences:default_path_current", "true", "If the current path is preferred over the default path.");
    defaults_.setValidStrings("preferences:default_path_current", ListUtils::create<String>("true,false"));
    defaults_.setValue("preferences:intensity_cutoff", "off", "Low intensity cutoff for maps.");
    defaults_.setValidStrings("preferences:intensity_cutoff", ListUtils::create<String>("on,off"));
    defaults_.setValue("preferences:on_file_change", "ask", "What action to take, when a data file changes. Do nothing, update automatically or ask the user.");
    defaults_.setValidStrings("preferences:on_file_change", ListUtils::create<String>("none,ask,update automatically"));
    defaults_.setValue("preferences:topp_cleanup", "true", "If the temporary files for calling of TOPP tools should be removed after the call.");
    defaults_.setValidStrings("preferences:topp_cleanup", ListUtils::create<String>("true,false"));
    defaults_.setValue("preferences:use_cached_ms2", "false", "If possible, only load MS1 spectra into memory and keep MS2 spectra on disk (using indexed mzML).");
    defaults_.setValidStrings("preferences:use_cached_ms2", ListUtils::create<String>("true,false"));
    defaults_.setValue("preferences:use_cached_ms1", "false", "If possible, do not load MS1 spectra into memory spectra into memory and keep MS2 spectra on disk (using indexed mzML).");
    defaults_.setValidStrings("preferences:use_cached_ms1", ListUtils::create<String>("true,false"));
    // 1d view
    defaults_.insert("preferences:1d:", Spectrum1DCanvas(Param()).getDefaults());
    defaults_.setSectionDescription("preferences:1d", "Settings for single spectrum view.");
    // 2d view
    defaults_.insert("preferences:2d:", Spectrum2DCanvas(Param()).getDefaults());
    defaults_.setSectionDescription("preferences:2d", "Settings for 2D map view.");
    // 3d view
    defaults_.insert("preferences:3d:", Spectrum3DCanvas(Param()).getDefaults());
    defaults_.setSectionDescription("preferences:3d", "Settings for 3D map view.");
    // identification view
    defaults_.insert("preferences:idview:", SpectraIdentificationViewWidget(Param()).getDefaults());
    defaults_.setSectionDescription("preferences:idview", "Settings for identification view.");
    defaults_.setValue("preferences:version", "none", "OpenMS version, used to check if the TOPPView.ini is up-to-date");
    subsections_.push_back("preferences:RecentFiles");
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
    dlg.setParam(param_);

    // --------------------------------------------------------------------
    // Execute dialog and update parameter object with user modified values
    if (dlg.exec())
    {
      param_ = dlg.getParam();
      savePreferences();
    }
  }

  void TOPPViewBase::addDataFile(const String& filename, bool show_options, bool add_to_recent, String caption, UInt window_id, Size spectrum_id)
  {
    setCursor(Qt::WaitCursor);
    RAIICleanup cl([&]() { setCursor(Qt::ArrowCursor); }); // revert to ArrowCursor on exit

    String abs_filename = File::absolutePath(filename);

    // check if the file exists
    if (!File::exists(abs_filename))
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Open file error", String("The file '") + abs_filename + "' does not exist!");
      return;
    }

    // determine file type
    FileHandler fh;
    FileTypes::Type file_type = fh.getType(abs_filename);
    if (file_type == FileTypes::UNKNOWN)
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Open file error", String("Could not determine file type of '") + abs_filename + "'!");
      return;
    }

    // abort if file type unsupported
    if (!supported_types.contains(file_type))
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Open file error", String("The type '") + FileTypes::typeToName(file_type) + "' is not supported!");
      return;
    }

    //try to load data and determine if it's 1D or 2D data

    // create shared pointer to main data types
    FeatureMapType* feature_map = new FeatureMapType();
    FeatureMapSharedPtrType feature_map_sptr(feature_map);

    ExperimentSharedPtrType peak_map_sptr(new ExperimentType());

    ConsensusMapType* consensus_map = new ConsensusMapType();
    ConsensusMapSharedPtrType consensus_map_sptr(consensus_map);

    vector<PeptideIdentification> peptides;

    LayerData::DataType data_type;

    ODExperimentSharedPtrType on_disc_peaks(new OnDiscMSExperiment);

    bool cache_ms2_on_disc = ((String)param_.getValue("preferences:use_cached_ms2") == "true");
    bool cache_ms1_on_disc = ((String)param_.getValue("preferences:use_cached_ms1") == "true");

    try
    {
      if (file_type == FileTypes::FEATUREXML)
      {
        FeatureXMLFile().load(abs_filename, *feature_map);
        data_type = LayerData::DT_FEATURE;
      }
      else if (file_type == FileTypes::CONSENSUSXML)
      {
        ConsensusXMLFile().load(abs_filename, *consensus_map);
        data_type = LayerData::DT_CONSENSUS;
      }
      else if (file_type == FileTypes::IDXML)
      {
        vector<ProteinIdentification> proteins; // not needed later
        IdXMLFile().load(abs_filename, proteins, peptides);
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
        data_type = LayerData::DT_IDENT;
      }
      else if (file_type == FileTypes::MZIDENTML)
      {
        vector<ProteinIdentification> proteins; // not needed later
        MzIdentMLFile().load(abs_filename, proteins, peptides);
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
        data_type = LayerData::DT_IDENT;
      }
      else
      {
        bool parsing_success = false;
        if (file_type == FileTypes::MZML)
        {
          // Load index only and check success (is it indexed?)
          MzMLFile f;
          Internal::IndexedMzMLHandler indexed_mzml_file_;
          indexed_mzml_file_.openFile(filename);
          if ( indexed_mzml_file_.getParsingSuccess() && cache_ms2_on_disc)
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

            // peak_map_sptr = boost::static_pointer_cast<ExperimentSharedPtrType>(on_disc_peaks->getMetaData());
            peak_map_sptr = on_disc_peaks->getMetaData();

            for (Size k = 0; k < indexed_mzml_file_.getNrSpectra() && !cache_ms1_on_disc; k++)
            {
              if ( peak_map_sptr->getSpectrum(k).getMSLevel() == 1)
              {
                peak_map_sptr->getSpectrum(k) = on_disc_peaks->getSpectrum(k);
              }
            }
            for (Size k = 0; k < indexed_mzml_file_.getNrChromatograms() && !cache_ms2_on_disc; k++)
            {
              peak_map_sptr->getChromatogram(k) = on_disc_peaks->getChromatogram(k);
            }

            // Load at least one spectrum into memory (TOPPView assumes that at least one spectrum is in memory)
            if (cache_ms1_on_disc && peak_map_sptr->getNrSpectra() > 0) peak_map_sptr->getSpectrum(0) = on_disc_peaks->getSpectrum(0);
          }
        }

        // Load all data into memory
        if (!parsing_success)
        {
          fh.loadExperiment(abs_filename, *peak_map_sptr, file_type, ProgressLogger::GUI);
        }
        OPENMS_LOG_INFO << "INFO: done loading all " << std::endl;

        // a mzML file may contain both, chromatogram and peak data
        // -> this is handled in SpectrumCanvas::addLayer
        data_type = LayerData::DT_CHROMATOGRAM;
        if (peak_map_sptr->containsScanOfLevel(1))
        {
          data_type = LayerData::DT_PEAK;
        }
      }
    }
    catch (Exception::BaseException& e)
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Error while loading file:", e.what());
      return;
    }

    // sort for mz and update ranges of newly loaded data
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

    addData(feature_map_sptr, 
      consensus_map_sptr, 
      peptides, 
      peak_map_sptr, 
      on_disc_peaks, 
      data_type, 
      false, 
      show_options, 
      true, 
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
  }

  void TOPPViewBase::addData(FeatureMapSharedPtrType feature_map,
                             ConsensusMapSharedPtrType consensus_map,
                             vector<PeptideIdentification>& peptides,
                             ExperimentSharedPtrType peak_map,
                             ODExperimentSharedPtrType on_disc_peak_map,
                             LayerData::DataType data_type,
                             bool show_as_1d,
                             bool show_options,
                             bool as_new_window,
                             const String& filename,
                             const String& caption,
                             UInt window_id,
                             Size spectrum_id)
  {
    // initialize flags with defaults from the parameters
    bool maps_as_2d = ((String)param_.getValue("preferences:default_map_view") == "2d");
    bool maps_as_1d = false;
    bool use_intensity_cutoff = ((String)param_.getValue("preferences:intensity_cutoff") == "on");
    bool is_dia_data = false;

    // feature, consensus feature and identifications can be merged
    bool mergeable = ((data_type == LayerData::DT_FEATURE) ||
                      (data_type == LayerData::DT_CONSENSUS) ||
                      (data_type == LayerData::DT_IDENT));

    // only one peak spectrum? disable 2D as default
    if (peak_map->size() == 1) { maps_as_2d = false; }

    // set the window where (new layer) data could be opened in
    // get EnhancedTabBarWidget with given id
    EnhancedTabBarWidgetInterface* tab_bar_target = ws_.getWidget(window_id);

    // cast to SpectrumWidget
    SpectrumWidget* target_window = dynamic_cast<SpectrumWidget*>(tab_bar_target);

    if (tab_bar_target == nullptr)
    {
      target_window = getActiveSpectrumWidget();
    }
    else
    {
      as_new_window = false;
    }

    // create dialog no matter if it is shown or not. It is used to determine the flags.
    TOPPViewOpenDialog dialog(caption, as_new_window, maps_as_2d, use_intensity_cutoff, this);

    //disable opening in new window when there is no active window or feature/ID data is to be opened, but the current window is a 3D window
    if (target_window == nullptr || (mergeable && dynamic_cast<Spectrum3DWidget*>(target_window) != nullptr))
    {
      dialog.disableLocation(true);
    }

    //disable 1d/2d/3d option for feature/consensus/identification maps
    if (mergeable)
    {
      dialog.disableDimension(true);
    }

    //disable cutoff for feature/consensus/identification maps
    if (mergeable)
    {
      dialog.disableCutoff(false);
    }

    //enable merge layers if a feature layer is opened and there are already features layers to merge it to
    if (mergeable && target_window != nullptr) //TODO merge
    {
      SpectrumCanvas* open_canvas = target_window->canvas();
      Map<Size, String> layers;
      for (Size i = 0; i < open_canvas->getLayerCount(); ++i)
      {
        if (data_type == open_canvas->getLayer(i).type)
        {
          layers[i] = open_canvas->getLayer(i).getName();
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
    if (as_new_window) //new window
    {
      if (maps_as_1d) // 2d in 1d window
      {
        target_window = new Spectrum1DWidget(getSpectrumParameters(1), &ws_);
      }
      else if (maps_as_2d || mergeable) //2d or features/IDs
      {
        target_window = new Spectrum2DWidget(getSpectrumParameters(2), &ws_);
      }
      else // 3d
      {
        target_window = new Spectrum3DWidget(getSpectrumParameters(3), &ws_);
      }
    }

    if (merge_layer == -1) //add layer to the window
    {
      if (data_type == LayerData::DT_FEATURE) //features
      {
        if (!target_window->canvas()->addLayer(feature_map, filename))
        {
          return;
        }
      }
      else if (data_type == LayerData::DT_CONSENSUS) //consensus features
      {
        if (!target_window->canvas()->addLayer(consensus_map, filename))
          return;
      }
      else if (data_type == LayerData::DT_IDENT)
      {
        if (!target_window->canvas()->addLayer(peptides, filename))
          return;
      }
      else //peaks
      {
        if (!target_window->canvas()->addLayer(peak_map, on_disc_peak_map, filename))
          return;

        //calculate noise
        if (use_intensity_cutoff)
        {
          double cutoff = estimateNoiseFromRandomScans(*(target_window->canvas()->getCurrentLayer().getPeakData()), 1, 10, 80);
          DataFilters filters;
          filters.add(DataFilters::DataFilter(DataFilters::INTENSITY, DataFilters::GREATER_EQUAL, cutoff));
          target_window->canvas()->setFilters(filters);
        }
        else // no mower, hide zeros if wanted
        {
          if (target_window->canvas()->getCurrentLayer().getPeakData()->hasZeroIntensities(1))
          {
            statusBar()->showMessage("Note: Data contains zero values.\nA filter will be added to hide these values.\nYou can reenable data points with zero intensity by removing the filter.");
            DataFilters filters;
            filters.add(DataFilters::DataFilter(DataFilters::INTENSITY, DataFilters::GREATER_EQUAL, 0.001));
            target_window->canvas()->setFilters(filters);
          }
        }

        Spectrum1DWidget* open_1d_window = dynamic_cast<Spectrum1DWidget*>(target_window);
        if (open_1d_window)
        {
          open_1d_window->canvas()->activateSpectrum(spectrum_id);
        }
      }
    }
    else //merge feature/ID data into feature layer
    {
      Spectrum2DCanvas* canvas = qobject_cast<Spectrum2DCanvas*>(target_window->canvas());
      if (data_type == LayerData::DT_CONSENSUS)
      {
        canvas->mergeIntoLayer(merge_layer, consensus_map);
      }
      else if (data_type == LayerData::DT_FEATURE)
      {
        canvas->mergeIntoLayer(merge_layer, feature_map);
      }
      else if (data_type == LayerData::DT_IDENT)
      {
        canvas->mergeIntoLayer(merge_layer, peptides);
      }
    }

    if (as_new_window)
    {
      showSpectrumWidgetInWindow(target_window, caption);
    }

    // enable spectra view tab (not required anymore since selection_view_.update() will decide automatically)
    //selection_view_->show(SpectraSelectionTabs::SPECTRA_IDX);
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
    auto* sw = dynamic_cast<SpectrumWidget*>(ws_.getWidget(id));
    if (!sw) return;
    sw->setFocus(); // triggers layerActivated...
  }

  void TOPPViewBase::closeTab()
  {
    ws_.activeSubWindow()->close();
  }

  void TOPPViewBase::editMetadata()
  {
    SpectrumCanvas* canvas = getActiveCanvas();

    // warn if hidden layer => wrong layer selected...
    if (!canvas->getCurrentLayer().visible)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "The current layer is not visible", "Have you selected the right layer for this action?");
    }

    //show editable meta data dialog
    canvas->showMetaData(true);
  }

  void TOPPViewBase::layerStatistics()
  {
    getActiveSpectrumWidget()->showStatistics();
  }

  void TOPPViewBase::showStatusMessage(string msg, OpenMS::UInt time)
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

  void TOPPViewBase::showCursorStatusInvert(double mz, double rt)
  {
    // swap rt vs mz (for vertical projection)
    showCursorStatus(rt, mz);
  }

  void TOPPViewBase::showCursorStatus(double mz, double rt)
  {
    message_label_->setText("");
    if (mz == -1)
    {
      mz_label_->setText("m/z: ");
    }
    else if (boost::math::isinf(mz) || boost::math::isnan(mz))
    {
      mz_label_->setText("m/z: n/a");
    }
    else
    {
      mz_label_->setText((String("m/z: ") + String::number(mz, 6).fillLeft(' ', 8)).toQString());
    }

    if (rt == -1)
    {
      rt_label_->setText("RT: ");
    }
    else if (boost::math::isinf(rt) || boost::math::isnan(rt))
    {
      rt_label_->setText("RT: n/a");
    }
    else
    {
      rt_label_->setText((String("RT: ") + String::number(rt, 1).fillLeft(' ', 8)).toQString());
    }
    statusBar()->update();
  }

  void TOPPViewBase::resetZoom()
  {
    SpectrumWidget* w = getActiveSpectrumWidget();
    if (w != nullptr)
    {
      w->canvas()->resetZoom();
    }
  }

  void TOPPViewBase::setIntensityMode(int index)
  {
    SpectrumWidget* w = getActiveSpectrumWidget();
    if (w)
    {
      intensity_button_group_->button(index)->setChecked(true);
      w->setIntensityMode((OpenMS::SpectrumCanvas::IntensityModes)index);
    }
  }

  void TOPPViewBase::setDrawMode1D(int index)
  {
    Spectrum1DWidget* w = getActive1DWidget();
    if (w)
    {
      w->canvas()->setDrawMode((OpenMS::Spectrum1DCanvas::DrawModes)index);
    }
  }

  void TOPPViewBase::changeLabel(QAction* action)
  {
    bool set = false;

    //label type is selected
    for (Size i = 0; i < LayerData::SIZE_OF_LABEL_TYPE; ++i)
    {
      if (action->text().toStdString() == LayerData::NamesOfLabelType[i])
      {
        getActive2DWidget()->canvas()->setLabel(LayerData::LabelType(i));
        set = true;
      }
    }

    //button is simply pressed
    if (!set)
    {
      if (getActive2DWidget()->canvas()->getCurrentLayer().label == LayerData::L_NONE)
      {
        getActive2DWidget()->canvas()->setLabel(LayerData::L_INDEX);
        dm_label_2d_->menu()->actions()[1]->setChecked(true);
      }
      else
      {
        getActive2DWidget()->canvas()->setLabel(LayerData::L_NONE);
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
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::F_UNASSIGNED, false);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_PEPTIDEMZ, false);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_LABELS, false);
    }
    else if (action->text().toStdString() == "Show by precursor m/z")
    {
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::F_UNASSIGNED, true);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_PEPTIDEMZ, false);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_LABELS, false);
    }
    else if (action->text().toStdString() == "Show by peptide mass")
    {
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::F_UNASSIGNED, true);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_PEPTIDEMZ, true);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_LABELS, false);
    }
    else if (action->text().toStdString() == "Show label meta data")
    {
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::F_UNASSIGNED, true);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_PEPTIDEMZ, false);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_LABELS, true);
    }
    else // button is simply pressed
    {
      bool previous = getActive2DWidget()->canvas()->getLayerFlag(LayerData::F_UNASSIGNED);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::F_UNASSIGNED,
                                                  !previous);
      if (previous) // now: don't show
      {
        dm_unassigned_2d_->menu()->actions()[0]->setChecked(true);
      }
      else // now: show by precursor
      {
        dm_unassigned_2d_->menu()->actions()[1]->setChecked(true);
      }
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_PEPTIDEMZ, false);
    }

    updateToolBar();
  }

  void TOPPViewBase::changeLayerFlag(bool on)
  {
    QAction* action = qobject_cast<QAction*>(sender());
    if (Spectrum2DWidget* win = getActive2DWidget())
    {
      //peaks
      if (action == dm_precursors_2d_)
      {
        win->canvas()->setLayerFlag(LayerData::P_PRECURSORS, on);
      }
      //features
      else if (action == dm_hulls_2d_)
      {
        win->canvas()->setLayerFlag(LayerData::F_HULLS, on);
      }
      else if (action == dm_hull_2d_)
      {
        win->canvas()->setLayerFlag(LayerData::F_HULL, on);
      }
      //consensus features
      else if (action == dm_elements_2d_)
      {
        win->canvas()->setLayerFlag(LayerData::C_ELEMENTS, on);
      }
      // identifications
      else if (action == dm_ident_2d_)
      {
        win->canvas()->setLayerFlag(LayerData::I_PEPTIDEMZ, on);
      }
    }
  }

  void TOPPViewBase::updateBarsAndMenus()
  {
    //Update filter bar, spectrum bar and layer bar
    layerActivated();
    updateMenu();
  }

  void TOPPViewBase::updateToolBar()
  {
    SpectrumWidget* w = getActiveSpectrumWidget();

    if (w)
    {
      //set intensity mode
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
    Spectrum1DWidget* w1 = getActive1DWidget();
    if (w1)
    {
      //draw mode
      draw_group_1d_->button(w1->canvas()->getDrawMode())->setChecked(true);

      //show/hide toolbars and buttons
      tool_bar_1d_->show();
      tool_bar_2d_peak_->hide();
      tool_bar_2d_feat_->hide();
      tool_bar_2d_cons_->hide();
      tool_bar_2d_ident_->hide();
    }

    // 2D
    Spectrum2DWidget* w2 = getActive2DWidget();
    if (w2)
    {
      tool_bar_1d_->hide();
      // check if there is a layer before requesting data from it
      if (w2->canvas()->getLayerCount() > 0)
      {
        //peak draw modes
        if (w2->canvas()->getCurrentLayer().type == LayerData::DT_PEAK)
        {
          dm_precursors_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::P_PRECURSORS));
          tool_bar_2d_peak_->show();
          tool_bar_2d_feat_->hide();
          tool_bar_2d_cons_->hide();
          tool_bar_2d_ident_->hide();
        }
        //feature draw modes
        else if (w2->canvas()->getCurrentLayer().type == LayerData::DT_FEATURE)
        {
          dm_hulls_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::F_HULLS));
          dm_hull_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::F_HULL));
          dm_unassigned_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::F_UNASSIGNED));
          dm_label_2d_->setChecked(w2->canvas()->getCurrentLayer().label != LayerData::L_NONE);
          tool_bar_2d_peak_->hide();
          tool_bar_2d_feat_->show();
          tool_bar_2d_cons_->hide();
          tool_bar_2d_ident_->hide();
        }
        //consensus feature draw modes
        else if (w2->canvas()->getCurrentLayer().type == LayerData::DT_CONSENSUS)
        {
          dm_elements_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::C_ELEMENTS));
          tool_bar_2d_peak_->hide();
          tool_bar_2d_feat_->hide();
          tool_bar_2d_cons_->show();
          tool_bar_2d_ident_->hide();
        }
        else if (w2->canvas()->getCurrentLayer().type == LayerData::DT_IDENT)
        {
          dm_ident_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::I_PEPTIDEMZ));
          tool_bar_2d_peak_->hide();
          tool_bar_2d_feat_->hide();
          tool_bar_2d_cons_->hide();
          tool_bar_2d_ident_->show();
        }
      }
    }

    // 3D
    Spectrum3DWidget* w3 = getActive3DWidget();
    if (w3)
    {
      //show/hide toolbars and buttons
      tool_bar_1d_->hide();
      tool_bar_2d_peak_->hide();
      tool_bar_2d_feat_->hide();
      tool_bar_2d_cons_->hide();
      tool_bar_2d_ident_->hide();
    }
  }

  void TOPPViewBase::updateLayerBar()
  {
    layers_view_->update(getActiveSpectrumWidget());
  }

  void TOPPViewBase::updateViewBar()
  {
    selection_view_->update();
  }

  void TOPPViewBase::updateMenu()
  {
    FS_TV fs;
    // is there a canvas?
    if (getActiveCanvas() != nullptr)
    {
      fs |= TV_STATUS::HAS_CANVAS;
      // is there a layer?
      if (getActiveCanvas()->getLayerCount() != 0)  fs |= TV_STATUS::HAS_LAYER;
    }
    // is this a 1D view
    if (getActive1DWidget() != nullptr) fs |= TV_STATUS::IS_1D_VIEW;
    // are we in 1D mirror mode
    if (getActive1DWidget() && getActive1DWidget()->canvas()->mirrorModeActive()) fs |= TV_STATUS::HAS_MIRROR_MODE;
    // is there a TOPP tool running
    if (topp_.process == nullptr)  fs |= TV_STATUS::TOPP_IDLE;
    
    menu_.update(fs);
  }

  void TOPPViewBase::updateFilterBar()
  {
    SpectrumCanvas* canvas = getActiveCanvas();
    if (canvas == nullptr)
      return;

    if (canvas->getLayerCount() == 0)
      return;
    
    filter_list_->set(getActiveCanvas()->getCurrentLayer().filters);
  }

  void TOPPViewBase::layerFilterVisibilityChange(bool on)
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

  void TOPPViewBase::layerZoomChanged() // todo rename zoomothers
  {
    if (!zoom_together_) return;

    QList<QMdiSubWindow *> windows = ws_.subWindowList();
    if (!windows.count()) return;

    SpectrumWidget* w = getActiveSpectrumWidget();
    DRange<2> new_visible_area = w->canvas()->getVisibleArea();
    // only zoom if other window is also (not) a chromatogram
    bool sender_is_chrom = w->canvas()->getCurrentLayer().type == LayerData::DT_CHROMATOGRAM ||
                           w->canvas()->getCurrentLayer().chromatogram_flag_set();

    // go through all windows, adjust the visible area where necessary
    for (int i = 0; i < int(windows.count()); ++i)
    {
      SpectrumWidget* specwidg = qobject_cast<SpectrumWidget*>(windows.at(i)->widget());
      if (!specwidg) continue;

      bool is_chrom = specwidg->canvas()->getCurrentLayer().type == LayerData::DT_CHROMATOGRAM ||
                      specwidg->canvas()->getCurrentLayer().chromatogram_flag_set();
      if (is_chrom != sender_is_chrom) continue;
      // not the same dimensionality (e.g. Spectrum1DCanvas vs. 2DCanvas)
      if (w->canvas()->getName() != specwidg->canvas()->getName()) continue;

      specwidg->canvas()->setVisibleArea(new_visible_area);
    }
  }

  void TOPPViewBase::layerDeactivated()
  {

  }

  void TOPPViewBase::showSpectrumWidgetInWindow(SpectrumWidget* sw, const String& caption)
  {
    ws_.addSubWindow(sw);
    connect(sw->canvas(), &SpectrumCanvas::preferencesChange, this, &TOPPViewBase::updateLayerBar);
    connect(sw->canvas(), &SpectrumCanvas::layerActivated, this, &TOPPViewBase::layerActivated);
    connect(sw->canvas(), &SpectrumCanvas::layerModficationChange, this, &TOPPViewBase::updateLayerBar);
    connect(sw->canvas(), &SpectrumCanvas::layerZoomChanged, this, &TOPPViewBase::layerZoomChanged);
    connect(sw, &SpectrumWidget::sendStatusMessage, this, &TOPPViewBase::showStatusMessage);
    connect(sw, &SpectrumWidget::sendCursorStatus, this, &TOPPViewBase::showCursorStatus);
    connect(sw, &SpectrumWidget::dropReceived, this, &TOPPViewBase::copyLayer);

    // 1D spectrum specific signals
    Spectrum1DWidget* sw1 = qobject_cast<Spectrum1DWidget*>(sw);
    if (sw1 != nullptr)
    {
      connect(sw1, &Spectrum1DWidget::showCurrentPeaksAs2D, this, &TOPPViewBase::showCurrentPeaksAs2D);
      connect(sw1, &Spectrum1DWidget::showCurrentPeaksAs3D, this, &TOPPViewBase::showCurrentPeaksAs3D);
      connect(sw1, &Spectrum1DWidget::showCurrentPeaksAsIonMobility, this, &TOPPViewBase::showCurrentPeaksAsIonMobility);
      connect(sw1, &Spectrum1DWidget::showCurrentPeaksAsDIA, this, &TOPPViewBase::showCurrentPeaksAsDIA);
    }

    // 2D spectrum specific signals
    Spectrum2DWidget* sw2 = qobject_cast<Spectrum2DWidget*>(sw);
    if (sw2 != nullptr)
    {
      connect(sw2->getHorizontalProjection(), &Spectrum2DWidget::sendCursorStatus, this, &TOPPViewBase::showCursorStatus);
      connect(sw2->getVerticalProjection(), &Spectrum2DWidget::sendCursorStatus, this, &TOPPViewBase::showCursorStatusInvert);
      connect(sw2, CONNECTCAST(Spectrum2DWidget, showSpectrumAs1D, (int)), selection_view_, CONNECTCAST(SpectraSelectionTabs, showSpectrumAs1D, (int)));
      connect(sw2, &Spectrum2DWidget::showCurrentPeaksAs3D , this, &TOPPViewBase::showCurrentPeaksAs3D);
    }

    // 3D spectrum specific signals
    Spectrum3DWidget* sw3 = qobject_cast<Spectrum3DWidget*>(sw);
    if (sw3 != nullptr)
    {
      connect(sw3, &Spectrum3DWidget::showCurrentPeaksAs2D,this, &TOPPViewBase::showCurrentPeaksAs2D);
    }

    sw->setWindowTitle(caption.toQString());
    sw->addToTabBar(&tab_bar_, caption, true);
    
    //show first window maximized (only visible windows are in the list)
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

  void TOPPViewBase::showGoToDialog()
  {
    SpectrumWidget* w = getActiveSpectrumWidget();
    if (w)
    {
      w->showGoToDialog();
    }
  }

  EnhancedWorkspace* TOPPViewBase::getWorkspace()
  {
    return &ws_;
  }

  SpectrumWidget* TOPPViewBase::getActiveSpectrumWidget() const
  {
    if (!ws_.currentSubWindow())
    {
      return nullptr;
    }
    return qobject_cast<SpectrumWidget*>(ws_.currentSubWindow()->widget());
  }

  SpectrumCanvas* TOPPViewBase::getActiveCanvas() const
  {
    SpectrumWidget* sw = getActiveSpectrumWidget();
    if (sw == nullptr)
    {
      return nullptr;
    }
    return sw->canvas();
  }

  Spectrum1DWidget* TOPPViewBase::getActive1DWidget() const
  {
    return qobject_cast<Spectrum1DWidget*>(getActiveSpectrumWidget());
  }

  Spectrum2DWidget* TOPPViewBase::getActive2DWidget() const
  {
    return qobject_cast<Spectrum2DWidget*>(getActiveSpectrumWidget());
  }

  Spectrum3DWidget* TOPPViewBase::getActive3DWidget() const
  {
    return qobject_cast<Spectrum3DWidget*>(getActiveSpectrumWidget());
  }

  void TOPPViewBase::loadPreferences(String filename)
  {
    // compose default ini file path
    String default_ini_file = String(QDir::homePath()) + "/.TOPPView.ini";

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

      //apply preferences if they are of the current TOPPView version
      if (!error && tmp.exists("preferences:version") &&
          tmp.getValue("preferences:version").toString() == VersionInfo::getVersion())
      {
        try
        {
          setParameters(tmp);
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

      //set parameters to defaults when something is fishy with the parameters file
      if (error)
      {
        //reset parameters (they will be stored again when TOPPView quits)
        setParameters(Param());

        cerr << "The TOPPView preferences files '" << filename << "' was ignored. It is no longer compatible with this TOPPView version and will be replaced." << endl;
      }
    }
    else if (filename != default_ini_file)
    {
      cerr << "Unable to load INI File: '" << filename << "'" << endl;
    }
    param_.setValue("PreferencesFile", filename);

    //set the recent files
    Param p = param_.copy("preferences:RecentFiles");
    QStringList rfiles;
    for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
    {
      QString filename = it->value.toQString();
      if (File::exists(filename))
        rfiles.append(filename);
    }
    recent_files_.set(rfiles);
  }

  void TOPPViewBase::savePreferences()
  {
    // replace recent files
    param_.removeAll("preferences:RecentFiles");
    const QStringList& rfiles = recent_files_.get();
    for (int i = 0; i < rfiles.size(); ++i)
    {
      param_.setValue("preferences:RecentFiles:" + String(i), rfiles[i]);
    }

    // set version
    param_.setValue("preferences:version", VersionInfo::getVersion());

    // save only the subsection that begins with "preferences:"
    try
    {
      ParamXMLFile().store(string(param_.getValue("PreferencesFile")), param_.copy("preferences:"));
    }
    catch (Exception::UnableToCreateFile& /*e*/)
    {
      cerr << "Unable to create INI File: '" << string(param_.getValue("PreferencesFile")) << "'" << endl;
    }
  }

  QStringList TOPPViewBase::chooseFilesDialog_(const String& path_overwrite)
  {
    // store active sub window
    QMdiSubWindow* old_active = ws_.currentSubWindow();
    RAIICleanup clean([&]() { ws_.setActiveSubWindow(old_active); });

    QString open_path = current_path_.toQString();
    if (path_overwrite != "")
    {
      open_path = path_overwrite.toQString();
    }
    // we use the QT file dialog instead of using QFileDialog::Names(...)
    // On Windows and Mac OS X, this static function will use the native file dialog and not a QFileDialog,
    // which prevents us from doing GUI testing on it.
    QFileDialog dialog(this, "Open file(s)", open_path, supported_types.toFileDialogFilter(FileTypes::Filter::BOTH, true).toQString());
    dialog.setFileMode(QFileDialog::ExistingFiles);
    if (dialog.exec())
    {
       return dialog.selectedFiles();
    }
    return QStringList();
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
    //warn if hidden layer => wrong layer selected...
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();
    if (!layer.visible)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "The current layer is not visible", "Have you selected the right layer for this action?");
    }

    //create and store unique file name prefix for files
    topp_.file_name = File::getTempDirectory() + "/TOPPView_" + File::getUniqueName();
    if (!File::writable(topp_.file_name + "_ini"))
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Cannot create temporary file", String("Cannot write to '") + topp_.file_name + "'_ini!");
      return;
    }
    ToolsDialog tools_dialog(this, topp_.file_name + "_ini", current_path_, layer.type);

    if (tools_dialog.exec() == QDialog::Accepted)
    {
      //Store tool name, input parameter and output parameter
      topp_.tool = tools_dialog.getTool();
      topp_.in = tools_dialog.getInput();
      topp_.out = tools_dialog.getOutput();
      topp_.visible_area_only = visible_area_only;
      //run the tool
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
    //warn if hidden layer => wrong layer selected...
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();
    if (!layer.visible)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "The current layer is not visible", "Have you selected the right layer for this action?");
    }

    //run the tool
    runTOPPTool_();
  }

  void TOPPViewBase::runTOPPTool_()
  {
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();


    //delete old input and output file
    File::remove(topp_.file_name + "_in");
    File::remove(topp_.file_name + "_out");

    //test if files are writable
    if (!File::writable(topp_.file_name + "_in"))
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Cannot create temporary file", String("Cannot write to '") + topp_.file_name + "_in'!");
      return;
    }
    if (!File::writable(topp_.file_name + "_out"))
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Cannot create temporary file", String("Cannot write to '") + topp_.file_name + "'_out!");
      return;
    }

    //Store data
    topp_.layer_name = layer.getName();
    topp_.window_id = getActiveSpectrumWidget()->getWindowId();
    topp_.spectrum_id = layer.getCurrentSpectrumIndex();
    if (layer.type == LayerData::DT_PEAK  && !(layer.chromatogram_flag_set()))
    {
      MzMLFile f;
      f.setLogType(ProgressLogger::GUI);
      if (topp_.visible_area_only)
      {
        ExperimentType exp;
        getActiveCanvas()->getVisiblePeakData(exp);
        f.store(topp_.file_name + "_in", exp);
      }
      else
      {
        f.store(topp_.file_name + "_in", *layer.getPeakData());
      }
    }
    else if (layer.type == LayerData::DT_CHROMATOGRAM || layer.chromatogram_flag_set())
    {
      MzMLFile f;
      // This means we have chromatogram data, either as DT_CHROMATOGRAM or as
      // DT_PEAK with the chromatogram flag set. To run the TOPPTool we need to
      // remove the flag and add the newly generated layer as spectrum data
      // (otherwise we run into problems with SpectraViewWidget::updateEntries
      // which assumes that all chromatogram data has chromatograms).
      getActiveCanvas()->getCurrentLayer().remove_chromatogram_flag(); // removing the flag is not constant
      //getActiveCanvas()->getCurrentLayer().getPeakData()->setMetaValue("chromatogram_passed_through_TOPP", "true");

      f.setLogType(ProgressLogger::GUI);
      if (topp_.visible_area_only)
      {
        ExperimentType exp;
        getActiveCanvas()->getVisiblePeakData(exp);
        f.store(topp_.file_name + "_in", exp);
      }
      else
      {
        f.store(topp_.file_name + "_in", *layer.getPeakData());
      }
    }
    else if (layer.type == LayerData::DT_FEATURE)
    {
      if (topp_.visible_area_only)
      {
        FeatureMapType map;
        getActiveCanvas()->getVisibleFeatureData(map);
        FeatureXMLFile().store(topp_.file_name + "_in", map);
      }
      else
      {
        FeatureXMLFile().store(topp_.file_name + "_in", *layer.getFeatureMap());
      }
    }
    else
    {
      if (topp_.visible_area_only)
      {
        ConsensusMapType map;
        getActiveCanvas()->getVisibleConsensusData(map);
        ConsensusXMLFile().store(topp_.file_name + "_in", map);
      }
      else
      {
        ConsensusXMLFile().store(topp_.file_name + "_in", *layer.getConsensusMap());
      }
    }

    // compose argument list
    QStringList args;
    args << "-ini"
         << (topp_.file_name + "_ini").toQString()
         << QString("-%1").arg(topp_.in.toQString())
         << (topp_.file_name + "_in").toQString()
         << "-no_progress";
    if (topp_.out != "")
    {
      args << QString("-%1").arg(topp_.out.toQString())
           << (topp_.file_name + "_out").toQString();
    }

    // start log and show it
    log_->appendNewHeader(LogWindow::LogState::NOTICE, QString("Starting '%1'").arg(topp_.tool.toQString()), ""); // tool + args.join(" "));

    // initialize process
    topp_.process = new QProcess();
    topp_.process->setProcessChannelMode(QProcess::MergedChannels);

    // connect slots
    connect(topp_.process, &QProcess::readyReadStandardOutput, this, &TOPPViewBase::updateProcessLog);
    connect(topp_.process, CONNECTCAST(QProcess, finished, (int, QProcess::ExitStatus)), this, &TOPPViewBase::finishTOPPToolExecution);
    QString tool_executable;
    try
    {
      // find correct location of TOPP tool
      tool_executable = File::findSiblingTOPPExecutable(topp_.tool).toQString();
    }
    catch (Exception::FileNotFound& /*ex*/)
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Could not locate executable!", QString("Finding executable of TOPP tool '%1' failed. Please check your TOPP/OpenMS installation. Workaround: Add the bin/ directory to your PATH").arg(topp_.tool.toQString()));
      return;
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

      // re-enable Apply TOPP tool menues
      delete topp_.process;
      topp_.process = nullptr;
      updateMenu();

      return;
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
    else if (topp_.out != "")
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, QString("'%1' finished successfully").arg(topp_.tool.toQString()),
                      QString("Execution time: %1 ms").arg(topp_.timer.elapsed()));
      if (!File::readable(topp_.file_name + "_out"))
      {
        log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Cannot read TOPP output", String("Cannot read '") + topp_.file_name + "_out'!");
      }
      else
      {
        addDataFile(topp_.file_name + "_out", true, false, topp_.layer_name + " (" + topp_.tool + ")", topp_.window_id, topp_.spectrum_id);
      }
    }

    //clean up
    delete topp_.process;
    topp_.process = nullptr;
    updateMenu();

    //clean up temporary files
    if (param_.getValue("preferences:topp_cleanup") == "true")
    {
      File::remove(topp_.file_name + "_ini");
      File::remove(topp_.file_name + "_in");
      File::remove(topp_.file_name + "_out");
    }
  }

  const LayerData* TOPPViewBase::getCurrentLayer() const
  {
    SpectrumCanvas* canvas = getActiveCanvas();
    if (canvas == nullptr)
    {
      return nullptr;
    }
    return &(canvas->getCurrentLayer());
  }

  void TOPPViewBase::toggleProjections()
  {
    Spectrum2DWidget* w = getActive2DWidget();
    if (w)
    {
      //update minimum size before
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

  bool TOPPViewBase::annotateMS1FromMassFingerprinting_(const FeatureMap& identifications)
  {
    LayerData& layer = getActiveCanvas()->getCurrentLayer();
    if (layer.type == LayerData::DT_PEAK)
    {
      IDMapper im;
      Param p = im.getParameters();
      p.setValue("rt_tolerance", 30.0);
      im.setParameters(p);
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "Note", "Mapping matches with 30 sec tolerance and no m/z limit to spectra...");
      im.annotate((*layer.getPeakDataMuteable()), identifications, true, true);
      return true;
    }
    return false;
  }

  void TOPPViewBase::annotateWithID()
  {
    LayerData& layer = getActiveCanvas()->getCurrentLayer();
    // warn if hidden layer => wrong layer selected...
    if (!layer.visible)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "The current layer is not visible", "Have you selected the right layer for this action? Aborting.");
      return;
    }

    // load id data
    QString fname = QFileDialog::getOpenFileName(this,
                                                "Select protein/AMT identification data",
                                                current_path_.toQString(),
                                                "idXML files (*.idXML); mzIdentML files (*.mzid,*.mzIdentML); featureXML files (*.featureXML); all files (*.*)");

    if (fname.isEmpty()) return;

    FileTypes::Type type = FileHandler::getType(fname);
    if (type == FileTypes::FEATUREXML)
    {
      FeatureMap fm;
      FeatureXMLFile().load(fname, fm);

      // last protein ID must be from AccurateMassSearch (it gets appended there)
      String engine = "no protein identification section found";
      bool ams_ok = false;
      if (fm.getProteinIdentifications().size() > 0)
      {
        engine = fm.getProteinIdentifications().back().getSearchEngine();
        if (engine == "AccurateMassSearch")
        {
          annotateMS1FromMassFingerprinting_(fm);
          selection_view_->setTabEnabled(SpectraSelectionTabs::IDENT_IDX, true);
          ams_ok = true;
        }
      }
      if (!ams_ok)
      {
        QMessageBox::warning(this, "Error", (String("FeatureXML is currently only supported for files generated by the AccurateMassSearch tool (got '") + engine +  "', expected 'AccurateMassSearch'.").toQString());
        return;
      }
    }
    else if (type == FileTypes::IDXML || type == FileTypes::MZIDENTML)
    {
      vector<PeptideIdentification> identifications;
      vector<ProteinIdentification> protein_identifications;
      try
      {
        String document_id;
        if (type == FileTypes::MZIDENTML) MzIdentMLFile().load(fname, protein_identifications, identifications);
        else IdXMLFile().load(fname, protein_identifications, identifications, document_id);
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::warning(this, "Error", QString("Loading of idXML/mzIdentML file failed! (") + e.what() + ")");
        return;
      }

      IDMapper mapper;
      if (layer.type == LayerData::DT_PEAK)
      {
        Param p = mapper.getDefaults();
        p.setValue("rt_tolerance", 0.1, "RT tolerance (in seconds) for the matching");
        p.setValue("mz_tolerance", 1.0, "m/z tolerance (in ppm or Da) for the matching");
        p.setValue("mz_measure", "Da", "unit of 'mz_tolerance' (ppm or Da)");
        mapper.setParameters(p);
        mapper.annotate(*layer.getPeakDataMuteable(), identifications, protein_identifications, true);
        selection_view_->show(SpectraSelectionTabs::IDENT_IDX);
      }
      else if (layer.type == LayerData::DT_FEATURE)
      {
        mapper.annotate(*layer.getFeatureMap(), identifications, protein_identifications);
      }
      else
      {
        mapper.annotate(*layer.getConsensusMap(), identifications, protein_identifications);
      }
    }
    else // file type other than mzIdentML, idXML or featureXML
    {
      QMessageBox::warning(this, "Error", QString("Unknown file type. No annotation performed."));
      return;
    }

    log_->appendNewHeader(LogWindow::LogState::NOTICE, "Done", "Annotation of spectra finished. Open identification view to see results!");
    updateViewBar(); // todo: may not be needed...
  }

  void TOPPViewBase::showSpectrumGenerationDialog()
  {
    TheoreticalSpectrumGenerationDialog spec_gen_dialog;
    if (spec_gen_dialog.exec())
    {
      String seq_string(spec_gen_dialog.getSequence());
      if (seq_string == "")
      {
        QMessageBox::warning(this, "Error", "You must enter a peptide sequence!");
        return;
      }
      AASequence aa_sequence;
      try
      {
        aa_sequence = AASequence::fromString(seq_string);
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::warning(this, "Error", QString("Spectrum generation failed! (") + e.what() + ")");
        return;
      }

      PeakSpectrum spectrum;
      Param p = spec_gen_dialog.getParam();
      Int charge = p.getValue("charge");

      p.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");

      // these two are true by default, initialize to false here and set to true in the loop below
      p.setValue("add_y_ions", "false", "Add peaks of y-ions to the spectrum");
      p.setValue("add_b_ions", "false", "Add peaks of b-ions to the spectrum");

      // for losses, isotopes, abundant_immonium_ions see getParam
      if (p.getValue("has_A").toBool()) // "A-ions"
      {
        p.setValue("add_a_ions", "true", "Add peaks of a-ions to the spectrum");
      }
      if (p.getValue("has_B").toBool()) // "B-ions"
      {
        p.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
      }
      if (p.getValue("has_C").toBool()) // "C-ions"
      {
        p.setValue("add_c_ions", "true", "Add peaks of c-ions to the spectrum");
      }
      if (p.getValue("has_X").toBool()) // "X-ions"
      {
        p.setValue("add_x_ions", "true", "Add peaks of x-ions to the spectrum");
      }
      if (p.getValue("has_Y").toBool()) // "Y-ions"
      {
        p.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
      }
      if (p.getValue("has_Z").toBool()) // "Z-ions"
      {
        p.setValue("add_z_ions", "true", "Add peaks of z-ions to the spectrum");
      }

      TheoreticalSpectrumGenerator generator;
      generator.setParameters(p);

      try
      {
        generator.getSpectrum(spectrum, aa_sequence, charge, charge);
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::warning(this, "Error", QString("Spectrum generation failed! (") + e.what() + "). Please report this to the developers (specify what input you used)!");
        return;
      }

      // set precursor information
      vector<Precursor> precursors;
      Precursor precursor;
      precursor.setMZ(aa_sequence.getMonoWeight());
      precursor.setCharge(charge);
      precursors.push_back(precursor);
      spectrum.setPrecursors(precursors);
      spectrum.setMSLevel(2);

      PeakMap new_exp;
      new_exp.addSpectrum(spectrum);
      ExperimentSharedPtrType new_exp_sptr(new PeakMap(new_exp));
      FeatureMapSharedPtrType f_dummy(new FeatureMapType());
      ConsensusMapSharedPtrType c_dummy(new ConsensusMapType());
      ODExperimentSharedPtrType od_dummy(new OnDiscMSExperiment());
      vector<PeptideIdentification> p_dummy;
      addData(f_dummy, c_dummy, p_dummy, new_exp_sptr, od_dummy, LayerData::DT_CHROMATOGRAM, false, true, true, "", seq_string + QString(" (theoretical)"));

      // ensure spectrum is drawn as sticks
      draw_group_1d_->button(Spectrum1DCanvas::DM_PEAKS)->setChecked(true);
      setDrawMode1D(Spectrum1DCanvas::DM_PEAKS);
    }
  }

  void TOPPViewBase::showSpectrumAlignmentDialog()
  {
    Spectrum1DWidget* active_1d_window = getActive1DWidget();
    if (!active_1d_window || !active_1d_window->canvas()->mirrorModeActive())
    {
      return;
    }
    Spectrum1DCanvas* cc = active_1d_window->canvas();

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
    LayerData& layer = getActiveCanvas()->getCurrentLayer();
    ExperimentSharedPtrType exp_sptr = layer.getPeakDataMuteable();
    ODExperimentSharedPtrType od_exp_sptr = layer.getOnDiscPeakData();

    //open new 2D widget
    Spectrum2DWidget* w = new Spectrum2DWidget(getSpectrumParameters(2), &ws_);

    //add data
    if (!w->canvas()->addLayer(exp_sptr, od_exp_sptr, layer.filename))
    {
      return;
    }

    String caption = layer.getName();
    // remove 3D suffix added when opening data in 3D mode (see below showCurrentPeaksAs3D())
    if (caption.hasSuffix(CAPTION_3D_SUFFIX_))
    {
      caption = caption.prefix(caption.rfind(CAPTION_3D_SUFFIX_));
    }
    w->canvas()->setLayerName(w->canvas()->getCurrentLayerIndex(), caption);
    showSpectrumWidgetInWindow(w, caption);
    updateMenu();
  }

  void TOPPViewBase::showCurrentPeaksAsIonMobility()
  {
    double IM_BINNING = 1e5;

    const LayerData& layer = getActiveCanvas()->getCurrentLayer();

    // Get current spectrum
    auto spidx = layer.getCurrentSpectrumIndex();
    MSSpectrum tmps = layer.getCurrentSpectrum();

    if (!tmps.containsIMData())
    {
      std::cout << "Cannot display ion mobility data, no float array with the correct name 'Ion Mobility' available." <<
        " Number of float arrays: " << tmps.getFloatDataArrays().size() << std::endl;
      return;
    }

    // Fill temporary spectral map (mobility -> Spectrum) with data from current spectrum
    std::map< int, boost::shared_ptr<MSSpectrum> > im_map;
    auto im_arr = tmps.getFloatDataArrays()[0]; // the first array should be the IM array (see containsIMData)
    for (Size k = 0;  k < tmps.size(); k++)
    {
      double im = im_arr[ k ];
      if (im_map.find( int(im*IM_BINNING) ) == im_map.end() )
      {
        boost::shared_ptr<MSSpectrum> news(new OpenMS::MSSpectrum() );
        news->setRT(im);
        news->setMSLevel(1);
        im_map[ int(im*IM_BINNING) ] = news;
      }
      im_map[ int(im*IM_BINNING) ]->push_back( tmps[k] );
    }

    // Add spectra into a MSExperiment, sort and prepare it for display
    ExperimentSharedPtrType tmpe(new OpenMS::MSExperiment() );
    for (const auto& s : im_map)
    {
      tmpe->addSpectrum( *(s.second) );
    }
    tmpe->sortSpectra();
    tmpe->updateRanges();
    tmpe->setMetaValue("is_ion_mobility", "true");
    tmpe->setMetaValue("ion_mobility_unit", "ms");

    // open new 2D widget
    Spectrum2DWidget* w = new Spectrum2DWidget(getSpectrumParameters(2), &ws_);

    // add data
    if (!w->canvas()->addLayer(tmpe, SpectrumCanvas::ODExperimentSharedPtrType(new OnDiscMSExperiment()), layer.filename))
    {
      return;
    }
    w->xAxis()->setLegend(SpectrumWidget::IM_MS_AXIS_TITLE);

    if (im_arr.getName().find("1002815") != std::string::npos)
    {
      w->xAxis()->setLegend(SpectrumWidget::IM_ONEKZERO_AXIS_TITLE);
      tmpe->setMetaValue("ion_mobility_unit", "1/K0");
    }

    String caption = layer.getName() + " (Ion Mobility Scan " + String(spidx) + ")";
    // remove 3D suffix added when opening data in 3D mode (see below showCurrentPeaksAs3D())
    if (caption.hasSuffix(CAPTION_3D_SUFFIX_))
    {
      caption = caption.prefix(caption.rfind(CAPTION_3D_SUFFIX_));
    }
    w->canvas()->setLayerName(w->canvas()->getCurrentLayerIndex(), caption);
    showSpectrumWidgetInWindow(w, caption);
    updateMenu();
  }

  void TOPPViewBase::showCurrentPeaksAsDIA()
  {
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();

    if (!layer.isDIAData())
    {
      std::cout << "Layer does not contain DIA / SWATH-MS data" << std::endl;
      return;
    }

    // Get current spectrum
    MSSpectrum tmps = layer.getCurrentSpectrum();

    // Add spectra into a MSExperiment, sort and prepare it for display
    ExperimentSharedPtrType tmpe(new OpenMS::MSExperiment() );

    // Collect all MS2 spectra with the same precursor as the current spectrum
    // (they are in the same SWATH window)
    String caption_add = "";
    if (!tmps.getPrecursors().empty())
    {
      // Get precursor isolation windows
      const auto& prec = tmps.getPrecursors()[0];
      double lower = prec.getMZ() - prec.getIsolationWindowLowerOffset();
      double upper = prec.getMZ() + prec.getIsolationWindowUpperOffset();

      Size k = 0;
      for (const auto& spec : (*layer.getPeakData() ) )
      {
        if (spec.getMSLevel() == 2 && !spec.getPrecursors().empty() )
        {
          if (fabs(spec.getPrecursors()[0].getMZ() - tmps.getPrecursors()[0].getMZ() ) < 1e-4 )
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
            else if (layer.getOnDiscPeakData()->getNrSpectra() > k)
            {
              // Get data from disk - copy data and tell TOPPView that this is
              // MS1 data so that it will be displayed properly in 2D and 3D
              // view
              MSSpectrum t = layer.getOnDiscPeakData()->getSpectrum(k);
              t.setMSLevel(1);
              tmpe->addSpectrum(t);
            }
          }
        }
        k++;
      }
      caption_add = "(DIA window " + String(lower) + " - " + String(upper) + ")";
    }
    
    tmpe->sortSpectra();
    tmpe->updateRanges();

    // open new 2D widget
    Spectrum2DWidget* w = new Spectrum2DWidget(getSpectrumParameters(2), &ws_);

    // add data
    if (!w->canvas()->addLayer(tmpe, SpectrumCanvas::ODExperimentSharedPtrType(new OnDiscMSExperiment()), layer.filename))
    {
      return;
    }

    String caption = layer.getName();
    caption += caption_add;
    // remove 3D suffix added when opening data in 3D mode (see below showCurrentPeaksAs3D())
    if (caption.hasSuffix(CAPTION_3D_SUFFIX_))
    {
      caption = caption.prefix(caption.rfind(CAPTION_3D_SUFFIX_));
    }
    w->canvas()->setLayerName(w->canvas()->getCurrentLayerIndex(), caption);
    showSpectrumWidgetInWindow(w, caption);
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
      if ((LayerData::DT_PEAK == getActiveCanvas()->getLayer(i).type) && // supported type
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

    LayerData& layer = const_cast<LayerData&>(getActiveCanvas()->getLayer(best_candidate));

    if (layer.type != LayerData::DT_PEAK)
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "Wrong layer type", "Something went wrong during layer selection. Please report this problem with a description of your current layers!");
    }
    //open new 3D widget
    Spectrum3DWidget* w = new Spectrum3DWidget(getSpectrumParameters(3), &ws_);

    ExperimentSharedPtrType exp_sptr = layer.getPeakDataMuteable();

    if (layer.isIonMobilityData())
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

    if (!w->canvas()->addLayer(exp_sptr, SpectrumCanvas::ODExperimentSharedPtrType(new OnDiscMSExperiment()), layer.filename))
    {
      return;
    }

    if (getActive1DWidget()) // switch from 1D to 3D
    {
      //TODO:
      //- doesnt make sense for fragment scan
      //- build new Area with mz range equal to 1D visible range
      //- rt range either overall MS1 data range or some convenient window

    }
    else if (getActive2DWidget()) // switch from 2D to 3D
    {
      w->canvas()->setVisibleArea(getActiveCanvas()->getVisibleArea());
    }

    // set layer name
    String caption = layer.getName() + CAPTION_3D_SUFFIX_;
    w->canvas()->setLayerName(w->canvas()->getCurrentLayerIndex(), caption);
    showSpectrumWidgetInWindow(w, caption);

    // set intensity mode (after spectrum has been added!)
    setIntensityMode(SpectrumCanvas::IM_SNAP);
    updateMenu();
  }

  void TOPPViewBase::updateProcessLog()
  {
    log_->appendText(topp_.process->readAllStandardOutput());
  }

  Param TOPPViewBase::getSpectrumParameters(UInt dim)
  {
    Param out = param_.copy(String("preferences:") + dim + "d:", true);
    out.setValue("default_path", param_.getValue("preferences:default_path").toString());
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
      topp_.process = nullptr;

      //finish log with new line
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
    for (StringList::const_iterator it = list.begin(); it != list.end(); ++it)
    {
      if (*it == "+")
      {
        last_was_plus = true;
        continue;
      }
      else if (std::find(colors.begin(), colors.end(), *it) != colors.end())
      {
        if ((getActive2DWidget() != nullptr || getActive3DWidget() != nullptr) && getActiveCanvas() != nullptr)
        {
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
          tmp.setValue("dot:gradient", gradients[Helpers::indexOf(colors, *it)]);
          getActiveCanvas()->setCurrentLayerParameters(tmp);
        }
      }
      else if (!last_was_plus || !getActiveSpectrumWidget())
      { // create new tab
        splash_screen->showMessage((String("Loading file: ") + *it).toQString());
        splash_screen->repaint();
        QApplication::processEvents();
        addDataFile(*it, false, true); // add data file but don't show options
      }
      else
      { // add to current tab
        splash_screen->showMessage((String("Loading file: ") + *it).toQString());
        splash_screen->repaint();
        QApplication::processEvents();
        last_was_plus = false;
        addDataFile(*it, false, true, "", getActiveSpectrumWidget()->getWindowId());
      }
    }
  }

  void TOPPViewBase::saveLayerAll()
  {
    getActiveCanvas()->saveCurrentLayer(false);
  }

  void TOPPViewBase::saveLayerVisible()
  {
    getActiveCanvas()->saveCurrentLayer(true);
  }

  void TOPPViewBase::toggleGridLines()
  {
    getActiveCanvas()->showGridLines(!getActiveCanvas()->gridLinesShown());
  }

  void TOPPViewBase::toggleAxisLegends()
  {
    getActiveSpectrumWidget()->showLegend(!getActiveSpectrumWidget()->isLegendShown());
  }

  void TOPPViewBase::toggleInterestingMZs()
  {
    auto w = getActive1DWidget();
    if (w == nullptr) return;
    w->canvas()->setDrawInterestingMZs(!w->canvas()->isDrawInterestingMZs());
  }

  void TOPPViewBase::showPreferences()
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
        if (!fh.loadExperiment(*it, exp))
        {
          QMessageBox::critical(this, "Error", "Only raw data files (mzML, DTA etc) are supported to view their meta data.");
          return;
        }
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


  void TOPPViewBase::showSpectrumMetaData(int spectrum_index)
  {
    getActiveCanvas()->showMetaData(true, spectrum_index);
  }

  void TOPPViewBase::copyLayer(const QMimeData* data, QWidget* source, int id)
  {
    SpectraViewWidget* spec_view = (source ? qobject_cast<SpectraViewWidget*>(source->parentWidget()) : nullptr);
    try
    {
      //NOT USED RIGHT NOW, BUT KEEP THIS CODE (it was hard to find out how this is done)
      //decode data to get the row
      //QByteArray encoded_data = data->data(data->formats()[0]);
      //QDataStream stream(&encoded_data, QIODevice::ReadOnly);
      //int row, col;
      //stream >> row >> col;

      // set wait cursor
      setCursor(Qt::WaitCursor);
      RAIICleanup cl([&]() { setCursor(Qt::ArrowCursor); });

      // determine where to copy the data
      UInt new_id = (id == -1) ? 0 : id;

      if (source == layers_view_)
      {
        // only the selected row can be dragged => the source layer is the selected layer
        LayerData& layer = getActiveCanvas()->getCurrentLayer();

        // attach feature, consensus and peak data
        FeatureMapSharedPtrType features = layer.getFeatureMap();
        ExperimentSharedPtrType peaks = layer.getPeakDataMuteable();
        ConsensusMapSharedPtrType consensus = layer.getConsensusMap();
        vector<PeptideIdentification> peptides = layer.peptides;
        ODExperimentSharedPtrType on_disc_peaks = layer.getOnDiscPeakData();

        // add the data
        addData(features, consensus, peptides, peaks, on_disc_peaks, layer.type, false, false, true, layer.filename, layer.getName(), new_id);
      }
      else if (spec_view != nullptr)
      {
        QTreeWidgetItem* item = spec_view->getTreeWidget()->currentItem();
        if (item != nullptr)
        {
          const LayerData& layer = getActiveCanvas()->getCurrentLayer();
          Size index = (Size)(item->text(3).toInt());       // todo: wtf...
          const ExperimentType::SpectrumType spectrum = (*layer.getPeakData())[index];
          ExperimentSharedPtrType new_exp_sptr(new ExperimentType());
          new_exp_sptr->addSpectrum(spectrum);
          ODExperimentSharedPtrType od_dummy(new OnDiscMSExperiment());
          FeatureMapSharedPtrType f_dummy(new FeatureMapType());
          ConsensusMapSharedPtrType c_dummy(new ConsensusMapType());
          vector<PeptideIdentification> p_dummy;
          addData(f_dummy, c_dummy, p_dummy, new_exp_sptr, od_dummy, LayerData::DT_CHROMATOGRAM, false, false, true, layer.filename, layer.getName(), new_id);
        }
      }
      else if (source == nullptr)
      {
        // drag source is external
        if (data->hasUrls())
        {
          QList<QUrl> urls = data->urls();
          for (QList<QUrl>::const_iterator it = urls.begin(); it != urls.end(); ++it)
          {
            addDataFile(it->toLocalFile(), false, true, "", new_id);
          }
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
    //do not update if the user disabled this feature.
    if (param_.getValue("preferences:default_path_current") != "true")
    {
      return;
    }

    //reset
    current_path_ = param_.getValue("preferences:default_path");

    //update if the current layer has a path associated
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
    std::vector<std::pair<const SpectrumWidget*, Size> > needs_update;
    for (const auto& mdi_window : ws_.subWindowList())
    {
      const SpectrumWidget* sw = qobject_cast<const SpectrumWidget*>(mdi_window);
      if (sw == nullptr) return;

      Size lc = sw->canvas()->getLayerCount();
      // determine if widget stores one or more layers for the given filename (->needs update)
      for (Size j = 0; j != lc; ++j)
      {
        if (sw->canvas()->getLayer(j).filename == filename)
        {
          needs_update.push_back(std::pair<const SpectrumWidget*, Size>(sw, j));
        }
      }
    }

    if (needs_update.empty()) // no layer references data of filename
    {
      watcher_->removeFile(filename); // remove watcher
      return;
    }
    
    //std::cout << "Number of Layers that need update: " << needs_update.size() << std::endl;
    pair<const SpectrumWidget*, Size>& slp = needs_update[0];
    const SpectrumWidget* sw = slp.first;
    Size layer_index = slp.second;

    bool user_wants_update = false;
    if ((String)(param_.getValue("preferences:on_file_change")) == "update automatically") //automatically update
    {
      user_wants_update = true;
    }
    else if ((String)(param_.getValue("preferences:on_file_change")) == "ask") //ask the user if the layer should be updated
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
    LayerData& layer = const_cast<LayerData&>(sw->canvas()->getLayer(layer_index));
    // reload data
    if (layer.type == LayerData::DT_PEAK) //peak data
    {
      try
      {
        FileHandler().loadExperiment(layer.filename, *layer.getPeakDataMuteable());
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::critical(this, "Error", (String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
        layer.getPeakDataMuteable()->clear(true);
      }
      layer.getPeakDataMuteable()->sortSpectra(true);
      layer.getPeakDataMuteable()->updateRanges(1);
    }
    else if (layer.type == LayerData::DT_FEATURE) //feature data
    {
      try
      {
        FileHandler().loadFeatures(layer.filename, *layer.getFeatureMap());
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::critical(this, "Error", (String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
        layer.getFeatureMap()->clear(true);
      }
      layer.getFeatureMap()->updateRanges();
    }
    else if (layer.type == LayerData::DT_CONSENSUS) //consensus feature data
    {
      try
      {
        ConsensusXMLFile().load(layer.filename, *layer.getConsensusMap());
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::critical(this, "Error", (String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
        layer.getConsensusMap()->clear(true);
      }
      layer.getConsensusMap()->updateRanges();
    }
    else if (layer.type == LayerData::DT_CHROMATOGRAM) //chromatogram
    {
      //TODO CHROM
      try
      {
        FileHandler().loadExperiment(layer.filename, *layer.getPeakDataMuteable());
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::critical(this, "Error", (String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
        layer.getPeakDataMuteable()->clear(true);
      }
      layer.getPeakDataMuteable()->sortChromatograms(true);
      layer.getPeakDataMuteable()->updateRanges(1);
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

} //namespace OpenMS
