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
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/VISUAL/SpectraViewWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum3DWidget.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewPrefDialog.h>
#include <OpenMS/VISUAL/DIALOGS/SpectrumAlignmentDialog.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>

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
#include <QtWidgets/QListWidget>
#include <QtWidgets/QListWidgetItem>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <QPainter>
#include <QtWidgets/QSplashScreen>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTextEdit>
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

  TOPPViewBase::TOPPViewBase(QWidget* parent) :
    QMainWindow(parent),
    DefaultParamHandler("TOPPViewBase"),
    ws_(this),
    tab_bar_(this),
    identificationview_behavior_(this), // controller for spectra and identification view
    spectraview_behavior_(this)         
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
    connect(&tab_bar_, &EnhancedTabBar::currentIdChanged, this, &TOPPViewBase::enhancedWorkspaceWindowChanged);
    connect(&tab_bar_, &EnhancedTabBar::aboutToCloseId, this, &TOPPViewBase::closeByTab);
    connect(&tab_bar_, &EnhancedTabBar::dropOnWidget, [this](const QMimeData* data, QWidget* source){ this->copyLayer(data, source); });
    connect(&tab_bar_, &EnhancedTabBar::dropOnTab, this, &TOPPViewBase::copyLayer);
    box_layout->addWidget(&tab_bar_);

    connect(&ws_, &EnhancedWorkspace::subWindowActivated, this, &TOPPViewBase::updateBarsAndMenus);
    connect(&ws_, &EnhancedWorkspace::dropReceived, this, &TOPPViewBase::copyLayer);
    box_layout->addWidget(&ws_);

    //################## MENUS #################
    // File menu
    QMenu* file = new QMenu("&File", this);
    menuBar()->addMenu(file);
    file->addAction("&Open file", this, [&]() { openFileDialog(); }, Qt::CTRL + Qt::Key_O);
    file->addAction("Open &example file", [&](){ openFileDialog(File::getOpenMSDataPath() + "/examples/"); }, Qt::CTRL + Qt::Key_E);
    file->addAction("&Close", this, &TOPPViewBase::closeFile, Qt::CTRL + Qt::Key_W);
    file->addSeparator();

    // Meta data
    file->addAction("&Show meta data (file)", this, &TOPPViewBase::metadataFileDialog);
    file->addSeparator();

    // Recent files
    QMenu* recent_menu = new QMenu("&Recent files", this);
    recent_actions_.resize(20);
    for (Size i = 0; i < 20; ++i)
    {
      recent_actions_[i] = recent_menu->addAction("", this, &TOPPViewBase::openRecentFile);
      recent_actions_[i]->setVisible(false);
    }
    file->addMenu(recent_menu);

    file->addSeparator();
    file->addAction("&Preferences", this, &TOPPViewBase::preferencesDialog);
    file->addAction("&Quit", qApp, SLOT(quit()));

    // Tools menu
    QMenu* tools = new QMenu("&Tools", this);
    menuBar()->addMenu(tools);
    tools->addAction("&Select data range", this, &TOPPViewBase::showGoToDialog, Qt::CTRL + Qt::Key_G);
    tools->addAction("&Edit meta data", this, &TOPPViewBase::editMetadata, Qt::CTRL + Qt::Key_M);
    tools->addAction("&Statistics", this, &TOPPViewBase::layerStatistics);
    tools->addSeparator();

    tools->addAction("Apply TOPP tool (whole layer)", this, &TOPPViewBase::showTOPPDialog, Qt::CTRL + Qt::Key_T)->setData(false);
    tools->addAction("Apply TOPP tool (visible layer data)", this, &TOPPViewBase::showTOPPDialog, Qt::CTRL + Qt::SHIFT + Qt::Key_T)->setData(true);
    tools->addAction("Rerun TOPP tool", this, &TOPPViewBase::rerunTOPPTool, Qt::Key_F4);
    tools->addSeparator();
    tools->addAction("&Annotate with identification", this, &TOPPViewBase::annotateWithID, Qt::CTRL + Qt::Key_I);
    tools->addAction("Align spectra", this, &TOPPViewBase::showSpectrumAlignmentDialog);
    tools->addAction("Generate theoretical spectrum", this, &TOPPViewBase::showSpectrumGenerationDialog);

    // Layer menu
    QMenu* layer = new QMenu("&Layer", this);
    menuBar()->addMenu(layer);
    layer->addAction("Save all data", this, &TOPPViewBase::saveLayerAll, Qt::CTRL + Qt::Key_S);
    layer->addAction("Save visible data", this, &TOPPViewBase::saveLayerVisible, Qt::CTRL + Qt::SHIFT + Qt::Key_S);
    layer->addSeparator();
    layer->addAction("Show/hide grid lines", this, &TOPPViewBase::toggleGridLines, Qt::CTRL + Qt::Key_R);
    layer->addAction("Show/hide axis legends", this, &TOPPViewBase::toggleAxisLegends, Qt::CTRL + Qt::Key_L);
    layer->addAction("Show/hide automated m/z annotations", this, &TOPPViewBase::toggleInterestingMZs);
    layer->addSeparator();
    layer->addAction("Preferences", this, &TOPPViewBase::showPreferences);

    // Windows menu
    QMenu* windows = new QMenu("&Windows", this);
    menuBar()->addMenu(windows);
    windows->addAction("&Cascade", &ws_, &EnhancedWorkspace::cascadeSubWindows);
    windows->addAction("&Tile automatic", &ws_, &EnhancedWorkspace::tileSubWindows);
    windows->addAction(QIcon(":/tile_vertical.png"), "Tile &vertical", &ws_, &EnhancedWorkspace::tileVertical);
    windows->addAction(QIcon(":/tile_horizontal.png"), "Tile &horizontal", &ws_, &EnhancedWorkspace::tileHorizontal);
    linkZoom_action_ = windows->addAction("Link &Zoom", this, &TOPPViewBase::linkZoom);
    windows->addSeparator();

    // Help menu
    QMenu* help = new QMenu("&Help", this);
    menuBar()->addMenu(help);
    help->addAction(QWhatsThis::createAction(help));
    help->addSeparator();
    help->addAction("OpenMS website", [&]() { GUIHelpers::openURL("http://www.OpenMS.de"); });
    help->addAction("Tutorials and documentation", [&]() { GUIHelpers::openURL("html/index.html"); }, Qt::Key_F1);

    help->addSeparator();
    help->addAction("&About", [&]() {QApplicationTOPP::showAboutDialog(this, "TOPPView");});

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
    layers_view_ = new QListWidget(layer_dock_widget_);
    layers_view_->setWhatsThis("Layer bar<BR><BR>Here the available layers are shown. Left-click on a layer to select it.<BR>Layers can be shown and hidden using the checkboxes in front of the name.<BR> Renaming and removing a layer is possible through the context menu.<BR>Dragging a layer to the tab bar copies the layer.<BR>Double-clicking a layer open its preferences.<BR>You can use the 'PageUp' and 'PageDown' buttons to change the selected layer.");

    layer_dock_widget_->setWidget(layers_view_);
    layers_view_->setContextMenuPolicy(Qt::CustomContextMenu);
    layers_view_->setDragEnabled(true);
    connect(layers_view_, &QListWidget::currentRowChanged, this, &TOPPViewBase::layerSelectionChange);
    connect(layers_view_, &QListWidget::customContextMenuRequested, this, &TOPPViewBase::layerContextMenu);
    connect(layers_view_, &QListWidget::itemChanged, this, &TOPPViewBase::layerVisibilityChange);
    connect(layers_view_, &QListWidget::itemDoubleClicked, this, &TOPPViewBase::layerEdit);

    windows->addAction(layer_dock_widget_->toggleViewAction());

    // Views dock widget
    views_dockwidget_ = new QDockWidget("Views", this);
    views_dockwidget_->setObjectName("views_dock_widget");
    addDockWidget(Qt::BottomDockWidgetArea, views_dockwidget_);
    views_tabwidget_ = new QTabWidget(views_dockwidget_);
    views_dockwidget_->setWidget(views_tabwidget_);

    // Hook-up controller and views for spectra inspection
    spectra_view_widget_ = new SpectraViewWidget();
    connect(spectra_view_widget_, &SpectraViewWidget::showSpectrumMetaData, this, &TOPPViewBase::showSpectrumMetaData);
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, showSpectrumAs1D, (int)),              this, CONNECTCAST(TOPPViewBase, showSpectrumAs1D, (int)));
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, showSpectrumAs1D, (std::vector<int>)), this, CONNECTCAST(TOPPViewBase, showSpectrumAs1D, (std::vector<int>)));
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, spectrumSelected, (int)),              &spectraview_behavior_, CONNECTCAST(TOPPViewSpectraViewBehavior, activate1DSpectrum, (int)));
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, spectrumSelected, (std::vector<int>)), &spectraview_behavior_, CONNECTCAST(TOPPViewSpectraViewBehavior, activate1DSpectrum, (const std::vector<int>&)));
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, spectrumDoubleClicked, (int)),              this, CONNECTCAST(TOPPViewBase, showSpectrumAs1D, (int)));
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, spectrumDoubleClicked, (std::vector<int>)), this, CONNECTCAST(TOPPViewBase, showSpectrumAs1D, (std::vector<int>)));

    // Hook-up controller and views for identification inspection
    spectra_identification_view_widget_ = new SpectraIdentificationViewWidget(Param());
    connect(spectra_identification_view_widget_, &SpectraIdentificationViewWidget::spectrumDeselected, &identificationview_behavior_, &TOPPViewIdentificationViewBehavior::deactivate1DSpectrum);
    connect(spectra_identification_view_widget_, &SpectraIdentificationViewWidget::showSpectrumAs1D, this, CONNECTCAST(TOPPViewBase, showSpectrumAs1D, (int)));
    connect(spectra_identification_view_widget_, &SpectraIdentificationViewWidget::spectrumSelected, 
            &identificationview_behavior_, CONNECTCAST(TOPPViewIdentificationViewBehavior, activate1DSpectrum, (int,int,int)));
    connect(spectra_identification_view_widget_, &SpectraIdentificationViewWidget::requestVisibleArea1D, &identificationview_behavior_, &TOPPViewIdentificationViewBehavior::setVisibleArea1D);

    views_tabwidget_->addTab(spectra_view_widget_, spectra_view_widget_->objectName());
    views_tabwidget_->addTab(spectra_identification_view_widget_, spectra_identification_view_widget_->objectName());
    views_tabwidget_->setTabEnabled(0, false);
    views_tabwidget_->setTabEnabled(1, false);

    // switch between different view tabs
    connect(views_tabwidget_, &QTabWidget::currentChanged, this, &TOPPViewBase::viewChanged);
    connect(views_tabwidget_, &QTabWidget::tabBarDoubleClicked, this, &TOPPViewBase::viewTabwidgetDoubleClicked);

    // add hide/show option to dock widget
    windows->addAction(views_dockwidget_->toggleViewAction());

    // filter dock widget
    filter_dock_widget_ = new QDockWidget("Data filters", this);
    filter_dock_widget_->setObjectName("filter_dock_widget");
    addDockWidget(Qt::BottomDockWidgetArea, filter_dock_widget_);
    filter_list_ = new FilterList(filter_dock_widget_);
    connect(filter_list_, &FilterList::filterChanged, [&](const DataFilters& filter) {
      getActiveCanvas()->setFilters(filter);
    });
    filter_dock_widget_->setWidget(filter_list_);
    windows->addAction(filter_dock_widget_->toggleViewAction());

    // log window
    QDockWidget* log_bar = new QDockWidget("Log", this);
    log_bar->setObjectName("log_bar");
    addDockWidget(Qt::BottomDockWidgetArea, log_bar);
    log_ = new QTextEdit(log_bar);
    log_->setReadOnly(true);
    log_->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(log_, &QTextEdit::customContextMenuRequested, this, &TOPPViewBase::logContextMenu);
    log_bar->setWidget(log_);
    windows->addAction(log_bar->toggleViewAction());

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
    defaults_.setValue("preferences:tmp_file_path", QDir::tempPath(), "Path where temporary files can be created.");
    defaults_.setValue("preferences:number_of_recent_files", 15, "Number of recent files in the main menu.");
    defaults_.setMinInt("preferences:number_of_recent_files", 5);
    defaults_.setMaxInt("preferences:number_of_recent_files", 20);
    defaults_.setValue("preferences:legend", "show", "Legend visibility");
    defaults_.setValidStrings("preferences:legend", ListUtils::create<String>("show,hide"));
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

  float TOPPViewBase::estimateNoiseFromRandomMS1Scans(const ExperimentType& exp, UInt n_scans)
  {
    if (!exp.containsScanOfLevel(1))
    {
      return 0.0;
    }

    float noise = 0.0;
    UInt count = 0;
    srand(time(nullptr));
    while (count < n_scans)
    {
      UInt scan = (UInt)((double)rand() / ((double)(RAND_MAX)+1.0f) * (double)(exp.size() - 1));

      if (scan < exp.size() && exp[scan].getMSLevel() == 1 && exp[scan].size() != 0)
      {
        vector<float> tmp;
        tmp.reserve(exp[scan].size());
        for (SpectrumType::ConstIterator it = exp[scan].begin()
             ; it != exp[scan].end()
             ; ++it)
        {
          tmp.push_back(it->getIntensity());
        }
        std::sort(tmp.begin(), tmp.end());
        noise += tmp[(UInt)ceil((float)(tmp.size() - 1) / 1.25f)];
        ++count;
      }
    }
    return noise / (double)n_scans;
  }
  
  // static
  bool TOPPViewBase::hasPeptideIdentifications(const ExperimentType& map)
  {
    for (Size i = 0; i != map.size(); ++i)
    {
      if (!map[i].getPeptideIdentifications().empty())
      {
        return true;
      }
    }
    return false;
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

  std::set<String> TOPPViewBase::getFilenamesOfOpenFiles_()
  {
    set<String> filename_set;
    // iterate over all windows
    QList<QMdiSubWindow *> wl = ws_.subWindowList();
    for (int i = 0; i != ws_.subWindowList().count(); ++i)
    {
      QWidget* w = wl[i];
      // iterate over all widgets
      const SpectrumWidget* sw = qobject_cast<const SpectrumWidget*>(w);
      if (sw != nullptr)
      {
        Size lc = sw->canvas()->getLayerCount();
        // iterate over all layers
        for (Size j = 0; j != lc; ++j)
        {
          filename_set.insert(sw->canvas()->getLayer(j).filename);
        }
      }
    }
    return filename_set;
  }

  void TOPPViewBase::addDataFile(const String& filename, bool show_options, bool add_to_recent, String caption, UInt window_id, Size spectrum_id)
  {
    setCursor(Qt::WaitCursor);
    RAIICleanup cl([&]() { setCursor(Qt::ArrowCursor); }); // revert to ArrowCursor on exit

    String abs_filename = File::absolutePath(filename);

    // check if the file exists
    if (!File::exists(abs_filename))
    {
      showLogMessage_(LS_ERROR, "Open file error", String("The file '") + abs_filename + "' does not exist!");
      return;
    }

    // determine file type
    FileHandler fh;
    FileTypes::Type file_type = fh.getType(abs_filename);
    if (file_type == FileTypes::UNKNOWN)
    {
      showLogMessage_(LS_ERROR, "Open file error", String("Could not determine file type of '") + abs_filename + "'!");
      return;
    }

    // abort if file type unsupported
    if (file_type == FileTypes::INI)
    {
      showLogMessage_(LS_ERROR, "Open file error", String("The type '") + FileTypes::typeToName(file_type) + "' is not supported!");
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
          showLogMessage_(LS_WARNING, "While loading file:", msg);
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
          showLogMessage_(LS_WARNING, "While loading file:", msg);
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
        FileTypes::Type type;
        type = FileHandler::getType(filename);
        bool parsing_success = false;
        if (type == FileTypes::MZML)
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
      showLogMessage_(LS_ERROR, "Error while loading file:", e.what());
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
    EnhancedTabBarWidgetInterface* tab_bar_target = window_(window_id);

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
          double cutoff = estimateNoiseFromRandomMS1Scans(*(target_window->canvas()->getCurrentLayer().getPeakData()));
          //create filter
          DataFilters::DataFilter filter;
          filter.field = DataFilters::INTENSITY;
          filter.op = DataFilters::GREATER_EQUAL;
          filter.value = cutoff;
          ///add filter
          DataFilters filters;
          filters.add(filter);
          target_window->canvas()->setFilters(filters);
        }
        else // no mower, hide zeros if wanted
        {
          if (target_window->canvas()->getCurrentLayer().getPeakData()->hasZeroIntensities(1))
          {
            // create filter
            DataFilters::DataFilter filter;
            filter.field = DataFilters::INTENSITY;
            filter.op = DataFilters::GREATER_EQUAL;
            filter.value = 0.001;
            statusBar()->showMessage("Note: Data contains zero values.\nA filter will be added to hide these values.\nYou can reenable data points with zero intensity by removing the filter.");
            // add filter
            DataFilters filters;
            filters.add(filter);
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

    // enable spectra view tab
    views_tabwidget_->setTabEnabled(0, true);
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
    // get/correct number of recent files
    UInt number_of_recent_files = UInt(param_.getValue("preferences:number_of_recent_files"));
    if (number_of_recent_files > 20)
    {
      number_of_recent_files = 20;
      param_.setValue("preferences:number_of_recent_files", 20);
    }

    for (Size i = 0; i < 20; ++i)
    {
      if (i < (UInt)(recent_files_.size()))
      {
        recent_actions_[i]->setText(recent_files_[(int)i]);
        recent_actions_[i]->setVisible(true);
      }
      else
      {
        recent_actions_[i]->setVisible(false);
      }
    }
  }

  void TOPPViewBase::openRecentFile()
  {
    QAction* action = qobject_cast<QAction*>(sender());
    if (action)
    {
      QString filename = action->text();
      addDataFile(filename, true, true);
    }
  }


  EnhancedTabBarWidgetInterface* TOPPViewBase::window_(int id) const
  {
    // return window with window_id == id
    QList<QMdiSubWindow *> windows = ws_.subWindowList();

    // return the actual widget
    for (int i = 0; i < windows.size(); ++i)
    {
      EnhancedTabBarWidgetInterface* w = dynamic_cast<EnhancedTabBarWidgetInterface*>(windows.at(i)->widget());
      if (w != 0 && w->getWindowId() == id) { return w; }
    }
    return nullptr;
  }

  void TOPPViewBase::closeByTab(int id)
  {
    QWidget* w = dynamic_cast<QWidget*>(window_(id));
    if (w)
    {
      QMdiSubWindow* parent = qobject_cast<QMdiSubWindow*>(w->parentWidget());
      parent->close();
      updateMenu();
    }
  }

  void TOPPViewBase::enhancedWorkspaceWindowChanged(int id)
  {
    QWidget* w = dynamic_cast<QWidget*>(window_(id));
    if (!w) return;

    w->setFocus();
    SpectrumWidget* sw = dynamic_cast<SpectrumWidget*>(w);
    if (!sw) return // SpectrumWidget

    views_tabwidget_->setTabEnabled(0, true);
    // check if there is a layer before requesting data from it
    if (sw->canvas()->getLayerCount() == 0) return;

    const ExperimentType& map = *sw->canvas()->getCurrentLayer().getPeakData();
    if (hasPeptideIdentifications(map))
    {
      views_tabwidget_->setTabEnabled(1, true);
      if (dynamic_cast<Spectrum2DWidget*>(w))
      {
        views_tabwidget_->setCurrentIndex(0); // switch to scan tab for 2D widget
      }
      // cppcheck produces a false positive warning here -> ignore
      // cppcheck-suppress multiCondition
      else if (dynamic_cast<Spectrum1DWidget*>(w))
      {
        views_tabwidget_->setCurrentIndex(1); // switch to identification tab for 1D widget
      }
    }
    else
    {
      views_tabwidget_->setTabEnabled(1, false);
      views_tabwidget_->setCurrentIndex(0); // stay on scan view tab
    }
  }

  void TOPPViewBase::closeFile()
  {
    ws_.activeSubWindow()->close();
    updateMenu();
  }

  void TOPPViewBase::editMetadata()
  {
    SpectrumCanvas* canvas = getActiveCanvas();

    // warn if hidden layer => wrong layer selected...
    if (!canvas->getCurrentLayer().visible)
    {
      showLogMessage_(LS_NOTICE, "The current layer is not visible", "Have you selected the right layer for this action?");
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
        showLogMessage_(LS_ERROR, OPENMS_PRETTY_FUNCTION, "Button for intensity mode does not exist");
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
    // reset
    layers_view_->clear();
    SpectrumCanvas* cc = getActiveCanvas();
    if (cc == nullptr) { return; }

    // determine if this is a 1D view (for text color)
    bool is_1d_view = (dynamic_cast<Spectrum1DCanvas*>(cc) != nullptr);

    layers_view_->blockSignals(true);
    RAIICleanup cl([&]() { layers_view_->blockSignals(false); });

    for (Size i = 0; i < cc->getLayerCount(); ++i)
    {
      const LayerData& layer = cc->getLayer(i);

      // add item
      QListWidgetItem* item = new QListWidgetItem(layers_view_);
      QString name = layer.getDecoratedName().toQString();
      
      item->setText(name);
      item->setToolTip(layer.filename.toQString());

      if (is_1d_view)
      { 
        if (cc->getLayerCount() > 1)
        {
          QPixmap icon(7, 7);
          icon.fill(QColor(layer.param.getValue("peak_color").toQString()));
          item->setIcon(icon);
        }
      }
      else
      {  // 2D/3D map view
        switch (layer.type)
        {
         case LayerData::DT_PEAK:
           item->setIcon(QIcon(":/peaks.png"));
  	 break;
         case LayerData::DT_FEATURE:
           item->setIcon(QIcon(":/convexhull.png"));
         break;
         case LayerData::DT_CONSENSUS:
           item->setIcon(QIcon(":/elements.png"));
         break;
         default:
         break;
        }
      }

      item->setCheckState(layer.visible ? Qt::Checked : Qt::Unchecked);
      
      // highlight active item
      if (i == cc->activeLayerIndex())
      {
        layers_view_->setCurrentItem(item);
      }
    }
    
  }

  void TOPPViewBase::updateViewBar()
  {
    SpectrumCanvas* cc = getActiveCanvas();
    int layer_row = layers_view_->currentRow();

    if (layer_row == -1 || cc == nullptr)
    {
      if (spectra_view_widget_)
      {
        spectra_view_widget_->getTreeWidget()->clear();
        spectra_view_widget_->getComboBox()->clear();
      }

      if (spectra_identification_view_widget_)
      {
        spectra_identification_view_widget_->setLayer(nullptr);
        // remove all entries
        QTableWidget* w = spectra_identification_view_widget_->getTableWidget();
        for (int i = w->rowCount() - 1; i >= 0; --i)
        {
          w->removeRow(i);
        }
        for (int i = w->columnCount() - 1; i >= 0; --i)
        {
          w->removeColumn(i);
        }
        w->clear();
        views_tabwidget_->setTabEnabled(1, false);
        views_tabwidget_->setTabEnabled(0, true);
      }
      return;
    }

    if (spectra_view_widget_->isVisible())
    {
      spectra_view_widget_->updateEntries(cc->getCurrentLayer());
    }

    if (spectra_identification_view_widget_->isVisible())
    {
      if (&cc->getCurrentLayer() != spectra_identification_view_widget_->getLayer())
      {
        spectra_identification_view_widget_->setLayer(&cc->getCurrentLayer());
      }
    }
  }

  void TOPPViewBase::viewChanged(int tab_index)
  {
    // set new behavior
    if (views_tabwidget_->tabText(tab_index) == spectra_view_widget_->objectName())
    {
      identificationview_behavior_.deactivateBehavior(); // finalize old behavior
      layer_dock_widget_->show();
      filter_dock_widget_->show();
      spectraview_behavior_.activateBehavior(); // initialize new behavior
    }
    else if (views_tabwidget_->tabText(tab_index) == spectra_identification_view_widget_->objectName())
    {
      spectraview_behavior_.deactivateBehavior();
      layer_dock_widget_->show();
      filter_dock_widget_->show();
      if (getActive2DWidget()) // currently 2D window is open
      {
        showSpectrumAs1D(0);
      }
      identificationview_behavior_.activateBehavior();
    }
    else
    {
      cerr << "Error: tab_index " << tab_index << endl;
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    updateViewBar();
  }

  void TOPPViewBase::viewTabwidgetDoubleClicked(int tab_index)
  {
    if (!getActiveSpectrumWidget()) return;

    // double click on disabled identification view
    // enables it and creates an empty identification structure
    if (views_tabwidget_->tabText(tab_index) == "Identification view"
      && !views_tabwidget_-> isTabEnabled(tab_index))
    {
      views_tabwidget_->setTabEnabled(1, true); // enable identification view
      views_tabwidget_->setCurrentIndex(1); // switch to identification view

      spectraview_behavior_.deactivateBehavior();
      layer_dock_widget_->show();
      filter_dock_widget_->show();
      if (getActive2DWidget()) // currently 2D window is open
      {
        showSpectrumAs1D(0);
      }
      identificationview_behavior_.activateBehavior();
    }

    updateViewBar();
  }

  void TOPPViewBase::layerSelectionChange(int i)
  {
    // after adding a layer i is -1. TODO: check if this is the correct behaviour
    if (i != -1)
    {
      getActiveCanvas()->activateLayer(i); // emits layerActivated 
    }
  }

  void TOPPViewBase::layerContextMenu(const QPoint& pos)
  {
    QListWidgetItem* item = layers_view_->itemAt(pos);
    if (!item) return;

    int layer = layers_view_->row(item);
    QMenu* context_menu = new QMenu(layers_view_);
    context_menu->addAction("Rename", [&]() {
      QString name = QInputDialog::getText(this, "Rename layer", "Name:", QLineEdit::Normal, getActiveCanvas()->getLayerName(layer).toQString());
      if (name != "")
      {
        getActiveCanvas()->setLayerName(layer, name);
      }});
    context_menu->addAction("Delete", [&]() {getActiveCanvas()->removeLayer(layer);});

    QAction* new_action = nullptr;
    if (getActiveCanvas()->getLayer(layer).flipped)
    {
      new_action = context_menu->addAction("Flip upwards (1D)", [&]() {
        getActive1DWidget()->canvas()->flipLayer(layer);
        bool b = getActive1DWidget()->canvas()->flippedLayersExist();
        getActive1DWidget()->canvas()->setMirrorModeActive(b);
      });
    }
    else
    {
      new_action = context_menu->addAction("Flip downwards (1D)", [&]() {
        getActive1DWidget()->canvas()->flipLayer(layer);
        getActive1DWidget()->canvas()->setMirrorModeActive(true);
      });
    }
    if (!getActive1DWidget())
    {
      new_action->setEnabled(false);
    }

    context_menu->addSeparator();
    context_menu->addAction("Preferences", [&]() {
      getActiveCanvas()->showCurrentLayerPreferences();
    });

    context_menu->exec(layers_view_->mapToGlobal(pos));
    
    // Update tab bar and window title
    if (getActiveCanvas()->getLayerCount() != 0)
    {
      tab_bar_.setTabText(tab_bar_.currentIndex(), getActiveCanvas()->getLayer(0).name.toQString());
      getActiveSpectrumWidget()->setWindowTitle(getActiveCanvas()->getLayer(0).name.toQString());
    }
    else
    {
      tab_bar_.setTabText(tab_bar_.currentIndex(), "empty");
      getActiveSpectrumWidget()->setWindowTitle("empty");
    }

    updateBarsAndMenus();
  }

  void TOPPViewBase::logContextMenu(const QPoint& pos)
  {
    QMenu context_menu;
    context_menu.addAction("Clear", [&]() {
      log_->clear();
    });
    context_menu.exec(log_->mapToGlobal(pos));
  }


  void TOPPViewBase::layerEdit(QListWidgetItem* /*item*/)
  {
    getActiveCanvas()->showCurrentLayerPreferences();
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
      getActiveCanvas()->changeLayerFilterState(getActiveCanvas()->activeLayerIndex(), on);
    }
  }

  void TOPPViewBase::layerVisibilityChange(QListWidgetItem* item)
  {
    int layer;
    bool visible;
    layer = layers_view_->row(item);
    visible = getActiveCanvas()->getLayer(layer).visible;

    if (item->checkState() == Qt::Unchecked && visible)
    {
      getActiveCanvas()->changeVisibility(layer, false);
    }
    else if (item->checkState() == Qt::Checked && !visible)
    {
      getActiveCanvas()->changeVisibility(layer, true);
    }
  }

  void TOPPViewBase::updateTabBar(QMdiSubWindow* w)
  {
    if (w)
    {
      EnhancedTabBarWidgetInterface* tbw = dynamic_cast<EnhancedTabBarWidgetInterface*>(w->widget());
      Int window_id = tbw->getWindowId();
      tab_bar_.setCurrentId(window_id);
    }
  }


  void TOPPViewBase::linkZoom()
  {
    zoom_together_ = !zoom_together_;
    if (!zoom_together_)
    {
      linkZoom_action_->setText("Link &Zoom");
    }
    else
    {
      linkZoom_action_->setText("Unlink &Zoom");
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

  void TOPPViewBase::layerZoomChanged()
  {
    QList<QMdiSubWindow *> windows = ws_.subWindowList();
    if (!windows.count())
      return;

    if (!zoom_together_)
      return;

    SpectrumWidget* w = getActiveSpectrumWidget();

    // figure out which dimension the active widget has: 2D (MSExperiment) or 1D (Iontrace)
    // and get the corresponding RT values.
    Spectrum1DWidget* sw1 = qobject_cast<Spectrum1DWidget*>(w);
    Spectrum2DWidget* sw2 = qobject_cast<Spectrum2DWidget*>(w);
    Spectrum3DWidget* sw3 = qobject_cast<Spectrum3DWidget*>(w);
    int widget_dimension = -1;
    if (sw1 != nullptr)
    {
      widget_dimension = 1;
    }
    else if (sw2 != nullptr)
    {
      widget_dimension = 2;
    }
    else if (sw3 != nullptr)
    {
      // dont link 3D
      widget_dimension = 3;
      return;
    }
    else
    {
      // Could not cast into any widget.
      return;
    }

    // check if the calling layer is a chromatogram:
    // - either its type is DT_CHROMATOGRAM
    // - or its peak data has a metavalue called "is_chromatogram" that is set to true
    if (getActiveCanvas()->getCurrentLayer().type == LayerData::DT_CHROMATOGRAM ||
        (getActiveCanvas()->getCurrentLayer().getPeakData()->size() > 0 &&
         getActiveCanvas()->getCurrentLayer().getPeakData()->metaValueExists("is_chromatogram") &&
         getActiveCanvas()->getCurrentLayer().getPeakData()->getMetaValue("is_chromatogram").toBool()
        ))
    {
      double minRT = -1, maxRT = -1;

      // Get the corresponding RT values depending on whether it is 2D (MSExperiment) or 1D (Iontrace).
      if (widget_dimension == 1)
      {
        minRT = sw1->canvas()->getVisibleArea().minX();
        maxRT = sw1->canvas()->getVisibleArea().maxX();
      }
      else if (widget_dimension == 2)
      {
        minRT = sw2->canvas()->getVisibleArea().minY();
        maxRT = sw2->canvas()->getVisibleArea().maxY();
      }

      // go through all windows, adjust the visible area where necessary
      for (int i = 0; i < int(windows.count()); ++i)
      {
        DRange<2> visible_area;

        QMdiSubWindow* window = windows.at(i);
        SpectrumWidget* specwidg = qobject_cast<SpectrumWidget*>(window->widget());

        // Skip if its not a SpectrumWidget, if it is not a chromatogram or if the dimensions don't match.
        if (!specwidg)
        {
          continue;
        }
        if (!(specwidg->canvas()->getCurrentLayer().type == LayerData::DT_CHROMATOGRAM) &&
            !(specwidg->canvas()->getCurrentLayer().getPeakData()->size() > 0 &&
              specwidg->canvas()->getCurrentLayer().getPeakData()->metaValueExists("is_chromatogram") &&
              specwidg->canvas()->getCurrentLayer().getPeakData()->getMetaValue("is_chromatogram").toBool()
              ))
        {
          continue;
        }
        if (!(widget_dimension == 1 && qobject_cast<Spectrum1DWidget*>(specwidg)) &&
            !(widget_dimension == 2 && qobject_cast<Spectrum2DWidget*>(specwidg)))
        {
          continue;
        }

        visible_area = specwidg->canvas()->getVisibleArea();

        // if we found a min/max RT, change all windows of 1 dimension
        if (minRT != -1 && maxRT != -1 && qobject_cast<Spectrum1DWidget*>(window->widget()))
        {
          visible_area.setMinX(minRT);
          visible_area.setMaxX(maxRT);
        }
        specwidg->canvas()->setVisibleArea(visible_area);
      }
    }
    else
    {
      DRange<2> new_visible_area = w->canvas()->getVisibleArea();
      // go through all windows, adjust the visible area where necessary
      for (int i = 0; i < int(windows.count()); ++i)
      {

        QMdiSubWindow* window = windows.at(i);
        SpectrumWidget* specwidg = qobject_cast<SpectrumWidget*>(window->widget());

        // Skip if its not a SpectrumWidget, if it is a chromatogram or if the dimensions don't match.
        if (!specwidg)
        {
          continue;
        }
        if ((specwidg->canvas()->getCurrentLayer().type == LayerData::DT_CHROMATOGRAM) ||
            (specwidg->canvas()->getCurrentLayer().getPeakData()->size() > 0 &&
             specwidg->canvas()->getCurrentLayer().getPeakData()->metaValueExists("is_chromatogram") &&
             specwidg->canvas()->getCurrentLayer().getPeakData()->getMetaValue("is_chromatogram").toBool()
            ))
        {
          continue;
        }
        if (!(widget_dimension == 1 && qobject_cast<Spectrum1DWidget*>(specwidg)) &&
            !(widget_dimension == 2 && qobject_cast<Spectrum2DWidget*>(specwidg)))
        {
          continue;
        }
        specwidg->canvas()->setVisibleArea(new_visible_area);
      }
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
      connect(sw2, CONNECTCAST(Spectrum2DWidget, showSpectrumAs1D, (int)), this, CONNECTCAST(TOPPViewBase, showSpectrumAs1D, (int)));
      connect(sw2, &Spectrum2DWidget::showCurrentPeaksAs3D , this, &TOPPViewBase::showCurrentPeaksAs3D);
    }

    // 3D spectrum specific signals
    Spectrum3DWidget* sw3 = qobject_cast<Spectrum3DWidget*>(sw);
    if (sw3 != nullptr)
    {
      connect(sw3, &Spectrum3DWidget::showCurrentPeaksAs2D,this, &TOPPViewBase::showCurrentPeaksAs2D);
    }

    sw->setWindowTitle(caption.toQString());

    //add tab with id
    static int window_counter = 4711;

    sw->setWindowId(window_counter++);

    tab_bar_.addTab(caption.toQString(), sw->getWindowId());

    //connect slots and signals for removing the widget from the bar, when it is closed
    //- through the menu entry
    //- through the tab bar
    //- through the MDI close button
    connect(sw, &SpectrumWidget::aboutToBeDestroyed, &tab_bar_, &EnhancedTabBar::removeId);

    tab_bar_.setCurrentId(sw->getWindowId());

    //show first window maximized (only visible windows are in the list)
    if (ws_.subWindowList().count() == 1)
    {
      sw->showMaximized();
    }
    else
    {
      sw->show();
    }
    enhancedWorkspaceWindowChanged(sw->getWindowId());
  }

  void TOPPViewBase::showGoToDialog()
  {
    SpectrumWidget* w = getActiveSpectrumWidget();
    if (w)
    {
      getActiveSpectrumWidget()->showGoToDialog();
    }
  }

  EnhancedWorkspace* TOPPViewBase::getWorkspace()
  {
    return &ws_;
  }

  SpectrumWidget* TOPPViewBase::getActiveSpectrumWidget() const
  {
    if (!ws_.activeSubWindow())
    {
      return nullptr;
    }
    return qobject_cast<SpectrumWidget*>(ws_.activeSubWindow()->widget());
  }

  SpectrumCanvas* TOPPViewBase::getActiveCanvas() const
  {
    if (ws_.currentSubWindow() == nullptr)
    {
      return nullptr;
    }
    SpectrumWidget* sw = qobject_cast<SpectrumWidget*>(ws_.currentSubWindow()->widget());
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
    if (p.size() != 0)
    {
      for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
      {
        QString filename = it->value.toQString();
        if (File::exists(filename))
          recent_files_.append(filename);
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
      param_.setValue("preferences:RecentFiles:" + String(i), recent_files_[i]);
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

  QStringList TOPPViewBase::getFileList_(const String& path_overwrite)
  {
    // store active sub window
    QMdiSubWindow* old_active = ws_.activeSubWindow();
    
    String filter_all = "readable files (*.mzML *.mzXML *.mzData *.featureXML *.consensusXML *.idXML *.dta *.dta2d fid *.bz2 *.gz);;";
    String filter_single = "mzML files (*.mzML);;mzXML files (*.mzXML);;mzData files (*.mzData);;feature map (*.featureXML);;consensus feature map (*.consensusXML);;peptide identifications (*.idXML);;XML files (*.xml);;XMass Analysis (fid);;dta files (*.dta);;dta2d files (*.dta2d);;bzipped files (*.bz2);;gzipped files (*.gz);;all files (*)";

    QString open_path = current_path_.toQString();
    if (path_overwrite != "")
    {
      open_path = path_overwrite.toQString();
    }
    // we use the QT file dialog instead of using QFileDialog::Names(...)
    // On Windows and Mac OS X, this static function will use the native file dialog and not a QFileDialog,
    // which prevents us from doing GUI testing on it.
    QFileDialog dialog(this, "Open file(s)", open_path, (filter_all + filter_single).toQString());
    dialog.setFileMode(QFileDialog::ExistingFiles);
    QStringList file_names;

    if (dialog.exec())
    {
      file_names = dialog.selectedFiles();
    }

    // restore active sub window
    ws_.setActiveSubWindow(old_active);
    
    return file_names;
  }

  void TOPPViewBase::openFileDialog(const String& dir)
  {
    for (const QString& filename : getFileList_(dir))
    {
      addDataFile(filename, true, true);
    }
  }

  void TOPPViewBase::rerunTOPPTool()
  {
    //warn if hidden layer => wrong layer selected...
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();
    if (!layer.visible)
    {
      showLogMessage_(LS_NOTICE, "The current layer is not visible", "Have you selected the right layer for this action?");
    }

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
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();
    if (!layer.visible)
    {
      showLogMessage_(LS_NOTICE, "The current layer is not visible", "Have you selected the right layer for this action?");
    }

    //create and store unique file name prefix for files
    topp_.file_name = param_.getValue("preferences:tmp_file_path").toString() + "/TOPPView_" + File::getUniqueName();
    if (!File::writable(topp_.file_name + "_ini"))
    {
      showLogMessage_(LS_ERROR, "Cannot create temporary file", String("Cannot write to '") + topp_.file_name + "'_ini!");
      return;
    }
    ToolsDialog tools_dialog(this, topp_.file_name + "_ini", current_path_, getCurrentLayer()->type);

    if (tools_dialog.exec() == QDialog::Accepted)
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
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();


    //delete old input and output file
    File::remove(topp_.file_name + "_in");
    File::remove(topp_.file_name + "_out");

    //test if files are writable
    if (!File::writable(topp_.file_name + "_in"))
    {
      showLogMessage_(LS_ERROR, "Cannot create temporary file", String("Cannot write to '") + topp_.file_name + "_in'!");
      return;
    }
    if (!File::writable(topp_.file_name + "_out"))
    {
      showLogMessage_(LS_ERROR, "Cannot create temporary file", String("Cannot write to '") + topp_.file_name + "'_out!");
      return;
    }

    //Store data
    topp_.layer_name = layer.name;
    topp_.window_id = getActiveSpectrumWidget()->getWindowId();
    topp_.spectrum_id = layer.getCurrentSpectrumIndex();
    if (layer.type == LayerData::DT_PEAK  && !(layer.chromatogram_flag_set()))
    {
      MzMLFile f;
      f.setLogType(ProgressLogger::GUI);
      if (topp_.visible)
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
      if (topp_.visible)
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
      if (topp_.visible)
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
      if (topp_.visible)
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
    showLogMessage_(LS_NOTICE, QString("Starting '%1'").arg(topp_.tool.toQString()), ""); // tool + args.join(" "));

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
      showLogMessage_(LS_ERROR, "Could not locate executable!", QString("Finding executable of TOPP tool '%1' failed. Please check your TOPP/OpenMS installation. Workaround: Add the bin/ directory to your PATH").arg(topp_.tool.toQString()));
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
      showLogMessage_(LS_ERROR, QString("Failed to execute '%1'").arg(topp_.tool.toQString()), QString("Execution of TOPP tool '%1' failed with error: %2").arg(topp_.tool.toQString(), topp_.process->errorString()));

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
    //finish with new line
    log_->append("");

    String tmp_dir = param_.getValue("preferences:tmp_file_path").toString();

    if (topp_.process->exitStatus() == QProcess::CrashExit)
    {
      showLogMessage_(LS_ERROR, QString("Execution of '%1' not successful!").arg(topp_.tool.toQString()),
                      QString("The tool crashed during execution. If you want to debug this crash, check the input files in '%1'"
                              " or enable 'debug' mode in the TOPP ini file.").arg(tmp_dir.toQString()));
    }
    else if (topp_.out != "")
    {
      showLogMessage_(LS_NOTICE, QString("'%1' finished successfully").arg(topp_.tool.toQString()),
                      QString("Execution time: %1 ms").arg(topp_.timer.elapsed()));
      if (!File::readable(topp_.file_name + "_out"))
      {
        showLogMessage_(LS_ERROR, "Cannot read TOPP output", String("Cannot read '") + topp_.file_name + "_out'!");
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

  void TOPPViewBase::loadFile(QString filename)
  {
    addDataFile(String(filename), true, false);
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
      showLogMessage_(LS_NOTICE, "Note", "Mapping matches with 30 sec tolerance and no m/z limit to spectra...");
      im.annotate((*layer.getPeakDataMuteable()), identifications, true, true);
    }
    else
    {
      return false;
    }
    return true;
  }

  void TOPPViewBase::annotateWithID()
  {
    LayerData& layer = getActiveCanvas()->getCurrentLayer();
    // warn if hidden layer => wrong layer selected...
    if (!layer.visible)
    {
      showLogMessage_(LS_NOTICE, "The current layer is not visible", "Have you selected the right layer for this action? Aborting.");
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
          views_tabwidget_->setTabEnabled(1, true); // enable identification view
          ams_ok = true;
        }
      }
      if (!ams_ok)
      {
        QMessageBox::warning(this, "Error", (String("FeatureXML is currently only supported for files generated by the AccurateMassSearch tool (got '") + engine +  "', expected 'AccurateMassSearch'.").toQString());
        return;
      }
    }
    else if (type == FileTypes::IDXML)
    {
      vector<PeptideIdentification> identifications;
      vector<ProteinIdentification> protein_identifications;

      try
      {
        String document_id;
        IdXMLFile().load(fname, protein_identifications, identifications, document_id);
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::warning(this, "Error", QString("Loading of idXML file failed! (") + e.what() + ")");
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
        views_tabwidget_->setTabEnabled(1, true); // enable identification view
        views_tabwidget_->setCurrentIndex(1); // switch to identification view
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
    else if (type == FileTypes::MZIDENTML)
    {
      vector<PeptideIdentification> identifications;
      vector<ProteinIdentification> protein_identifications;

      try
      {
        MzIdentMLFile().load(fname, protein_identifications, identifications);
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::warning(this, "Error", QString("Loading of idXML file failed! (") + e.what() + ")");
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
        views_tabwidget_->setTabEnabled(1, true); // enable identification view
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

    showLogMessage_(LS_NOTICE, "Done", "Annotation of spectra finished. Open identification view to see results!");
    updateViewBar();
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

  void TOPPViewBase::showSpectrumAs1D(int index)
  {
    Spectrum1DWidget* widget_1d = getActive1DWidget();
    Spectrum2DWidget* widget_2d = getActive2DWidget();

    if (widget_1d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_behavior_.showSpectrumAs1D(index);
      }

      if (spectra_identification_view_widget_->isVisible())
      {
        identificationview_behavior_.showSpectrumAs1D(index);
      }
    }
    else if (widget_2d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_behavior_.showSpectrumAs1D(index);
      }

      if (spectra_identification_view_widget_->isVisible())
      {
        identificationview_behavior_.showSpectrumAs1D(index);
      }
    }
  }

  void TOPPViewBase::showSpectrumAs1D(std::vector<int, std::allocator<int> > indices)
  {
    Spectrum1DWidget* widget_1d = getActive1DWidget();
    Spectrum2DWidget* widget_2d = getActive2DWidget();

    if (widget_1d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_behavior_.showSpectrumAs1D(indices);
      }
    }
    else if (widget_2d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_behavior_.showSpectrumAs1D(indices);
      }
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

    String caption = layer.name;
    // remove 3D suffix added when opening data in 3D mode (see below showCurrentPeaksAs3D())
    if (caption.hasSuffix(CAPTION_3D_SUFFIX_))
    {
      caption = caption.prefix(caption.rfind(CAPTION_3D_SUFFIX_));
    }
    w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);
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

    String caption = layer.name + " (Ion Mobility Scan " + String(spidx) + ")";
    // remove 3D suffix added when opening data in 3D mode (see below showCurrentPeaksAs3D())
    if (caption.hasSuffix(CAPTION_3D_SUFFIX_))
    {
      caption = caption.prefix(caption.rfind(CAPTION_3D_SUFFIX_));
    }
    w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);
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

    String caption = layer.name;
    caption += caption_add;
    // remove 3D suffix added when opening data in 3D mode (see below showCurrentPeaksAs3D())
    if (caption.hasSuffix(CAPTION_3D_SUFFIX_))
    {
      caption = caption.prefix(caption.rfind(CAPTION_3D_SUFFIX_));
    }
    w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);
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
      showLogMessage_(LS_NOTICE, "No compatible layer",
          "No layer found which is supported by the 3D view.");
      return;
    }


    if (best_candidate != target_layer)
    {
      showLogMessage_(LS_NOTICE, "Auto-selected compatible layer",
          "The currently active layer cannot be viewed in 3D view. The closest layer which is supported by the 3D view was selected!");
    }

    LayerData& layer = const_cast<LayerData&>(getActiveCanvas()->getLayer(best_candidate));

    if (layer.type != LayerData::DT_PEAK)
    {
      showLogMessage_(LS_NOTICE, "Wrong layer type", "Something went wrong during layer selection. Please report this problem with a description of your current layers!");
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
    String caption = layer.name + CAPTION_3D_SUFFIX_;
    w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);
    showSpectrumWidgetInWindow(w, caption);

    // set intensity mode (after spectrum has been added!)
    setIntensityMode(SpectrumCanvas::IM_SNAP);
    updateMenu();
  }

  void TOPPViewBase::updateProcessLog()
  {
    //show log if there is output
    qobject_cast<QWidget*>(log_->parent())->show();

    //update log_
    log_->moveCursor(QTextCursor::End, QTextCursor::MoveAnchor); // move cursor to end, since text is inserted at cursor
    log_->insertPlainText(topp_.process->readAllStandardOutput());

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
      log_->append("");

      updateMenu();
    }
  }

  void TOPPViewBase::updateMenu()
  {
    //is there a canvas?
    bool canvas_exists = false;
    if (getActiveCanvas() != nullptr)
    {
      canvas_exists = true;
    }
    //is there a layer?
    bool layer_exists = false;
    if (canvas_exists && getActiveCanvas()->getLayerCount() != 0)
    {
      layer_exists = true;
    }
    //is there a TOPP tool running
    bool topp_running = false;
    if (topp_.process != nullptr)
    {
      topp_running = true;
    }

    bool mirror_mode = getActive1DWidget() && getActive1DWidget()->canvas()->mirrorModeActive();
    QList<QAction*> actions = this->findChildren<QAction*>("");
    for (int i = 0; i < actions.count(); ++i)
    {
      QString text = actions[i]->text();
      if (text == "&Close" || text == "Show/hide grid lines" || text == "Show/hide axis legends")
      {
        actions[i]->setEnabled(false);
        if (canvas_exists)
        {
          actions[i]->setEnabled(true);
        }
      }
      else if (text.left(15) == "Apply TOPP tool")
      {
        actions[i]->setEnabled(false);
        if (canvas_exists && layer_exists && !topp_running)
        {
          actions[i]->setEnabled(true);
        }
      }
      else if (text == "Abort running TOPP tool")
      {
        actions[i]->setEnabled(false);
        if (topp_running)
        {
          actions[i]->setEnabled(true);
        }
      }
      else if (text == "Rerun TOPP tool")
      {
        actions[i]->setEnabled(false);
        if (canvas_exists && layer_exists && !topp_running && topp_.tool != "")
        {
          actions[i]->setEnabled(true);
        }
      }
      else if (text == "&Select data range" || text == "&Edit meta data" || text == "&Statistics" || text == "&Annotate with identification"  || text == "Save all data"  || text == "Save visible data"  || text == "Preferences")
      {
        actions[i]->setEnabled(false);
        if (canvas_exists && layer_exists)
        {
          actions[i]->setEnabled(true);
        }
      }
      else if (text == "Align spectra")
      {
        actions[i]->setEnabled(false);
        if (mirror_mode)
        {
          actions[i]->setEnabled(true);
        }
      }
      else if (text == "")
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
    for (StringList::const_iterator it = list.begin(); it != list.end(); ++it)
    {
      if (*it == "+")
      {
        last_was_plus = true;
        continue;
      }
      else if (*it == "@bw")
      {
        if ((getActive2DWidget() != nullptr || getActive3DWidget() != nullptr) && getActiveCanvas() != nullptr)
        {
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
          tmp.setValue("dot:gradient", "Linear|0,#ffffff;100,#000000");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
        }
      }
      else if (*it == "@bg")
      {
        if ((getActive2DWidget() != nullptr || getActive3DWidget() != nullptr) && getActiveCanvas() != nullptr)
        {
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
          tmp.setValue("dot:gradient", "Linear|0,#dddddd;100,#000000");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
        }
      }
      else if (*it == "@b")
      {
        if ((getActive2DWidget() != nullptr || getActive3DWidget() != nullptr) && getActiveCanvas() != nullptr)
        {
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
          tmp.setValue("dot:gradient", "Linear|0,#000000;100,#000000");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
        }
      }
      else if (*it == "@r")
      {
        if ((getActive2DWidget() != nullptr || getActive3DWidget() != nullptr) && getActiveCanvas() != nullptr)
        {
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
          tmp.setValue("dot:gradient", "Linear|0,#ff0000;100,#ff0000");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
        }
      }
      else if (*it == "@g")
      {
        if ((getActive2DWidget() != nullptr || getActive3DWidget() != nullptr) && getActiveCanvas() != nullptr)
        {
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
          tmp.setValue("dot:gradient", "Linear|0,#00ff00;100,#00ff00");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
        }
      }
      else if (*it == "@m")
      {
        if ((getActive2DWidget() != nullptr || getActive3DWidget() != nullptr) && getActiveCanvas() != nullptr)
        {
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
          tmp.setValue("dot:gradient", "Linear|0,#ff00ff;100,#ff00ff");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
        }
      }
      else if (!last_was_plus || !getActiveSpectrumWidget())
      {
        splash_screen->showMessage((String("Loading file: ") + *it).toQString());
        splash_screen->repaint();
        QApplication::processEvents();
        addDataFile(*it, false, true); // add data file but don't show options
      }
      else
      {
        splash_screen->showMessage((String("Loading file: ") + *it).toQString());
        splash_screen->repaint();
        QApplication::processEvents();
        last_was_plus = false;
        addDataFile(*it, false, true, "", getActiveSpectrumWidget()->getWindowId());
      }
    }
  }

  void TOPPViewBase::showLogMessage_(TOPPViewBase::LogState state, const String& heading, const String& body)
  {
    //Compose current time string
    DateTime d = DateTime::now();

    String state_string;
    switch (state)
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
    qobject_cast<QWidget*>(log_->parent())->show();
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
    QStringList files = getFileList_();
    FileHandler fh;
    fh.getOptions().setMetadataOnly(true);
    for (QStringList::iterator it = files.begin(); it != files.end(); ++it)
    {
      ExperimentType exp;
      try
      {
        fh.loadExperiment(*it, exp);
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

  SpectraIdentificationViewWidget* TOPPViewBase::getSpectraIdentificationViewWidget()
  {
    return spectra_identification_view_widget_;
  }

  void TOPPViewBase::showSpectrumMetaData(int spectrum_index)
  {
    getActiveCanvas()->showMetaData(true, spectrum_index);
  }

  void TOPPViewBase::copyLayer(const QMimeData* data, QWidget* source, int id)
  {
    QTreeWidget* spectra_view_treewidget = spectra_view_widget_->getTreeWidget();
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
        addData(features, consensus, peptides, peaks, on_disc_peaks, layer.type, false, false, true, layer.filename, layer.name, new_id);
      }
      else if (source == spectra_view_treewidget)
      {
        const LayerData& layer = getActiveCanvas()->getCurrentLayer();
        QTreeWidgetItem* item = spectra_view_treewidget->currentItem();
        if (item != nullptr)
        {
          Size index = (Size)(item->text(3).toInt());
          const ExperimentType::SpectrumType spectrum = (*layer.getPeakData())[index];
          ExperimentType new_exp;
          new_exp.addSpectrum(spectrum);
          ExperimentSharedPtrType new_exp_sptr(new ExperimentType(new_exp));
          ODExperimentSharedPtrType od_dummy(new OnDiscMSExperiment());
          FeatureMapSharedPtrType f_dummy(new FeatureMapType());
          ConsensusMapSharedPtrType c_dummy(new ConsensusMapType());
          vector<PeptideIdentification> p_dummy;
          addData(f_dummy, c_dummy, p_dummy, new_exp_sptr, od_dummy, LayerData::DT_CHROMATOGRAM, false, false, true, layer.filename, layer.name, new_id);
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
      showLogMessage_(LS_ERROR, "Error while creating layer", e.what());
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

    QList<QMdiSubWindow *> wl = ws_.subWindowList();

    // iterate over all windows and determine which need an update
    std::vector<std::pair<const SpectrumWidget*, Size> > needs_update;
    for (int i = 0; i != ws_.subWindowList().count(); ++i)
    {
      //std::cout << "Number of windows: " << ws_.subWindowList().count() << std::endl;
      QWidget* w = wl[i];
      const SpectrumWidget* sw = qobject_cast<const SpectrumWidget*>(w);
      if (sw != nullptr)
      {
        Size lc = sw->canvas()->getLayerCount();

        // determine if widget stores one or more layers for the given filename (->needs update)
        for (Size j = 0; j != lc; ++j)
        {
          //std::cout << "Layer filename: " << sw->canvas()->getLayer(j).filename << std::endl;
          const LayerData& ld = sw->canvas()->getLayer(j);
          if (ld.filename == filename)
          {
            needs_update.push_back(std::pair<const SpectrumWidget*, Size>(sw, j));
          }
        }
      }
    }

    if (needs_update.empty()) // no layer references data of filename
    {
      watcher_->removeFile(filename); // remove watcher
      return;
    }
    else if (!needs_update.empty()) // at least one layer references data of filename
    {
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
      else //if (user_wants_update == true)
      {
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
      }

      // update all layers that need an update
      for (Size i = 0; i != needs_update.size(); ++i)
      {
        pair<const SpectrumWidget*, Size>& slp = needs_update[i];
        const SpectrumWidget* sw = slp.first;
        Size layer_index = slp.second;
        sw->canvas()->updateLayer(layer_index);
      }
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
