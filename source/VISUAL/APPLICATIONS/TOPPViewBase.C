// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DB/DBAdapter.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/VISUAL/DIALOGS/DataFilterDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/DBOpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewPrefDialog.h>
#include <OpenMS/VISUAL/DIALOGS/SpectrumAlignmentDialog.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/SpectraViewWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum3DWidget.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>
#include <OpenMS/VISUAL/EnhancedWorkspace.h>

// OpenMS TOPPAS
#include <OpenMS/VISUAL/APPLICATIONS/TOPPASBase.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>


//Qt
#include <QtCore/QDate>
#include <QtCore/QDir>
#include <QtCore/QTime>
#include <QtCore/QUrl>
#include <QtGui/QCheckBox>
#include <QtGui/QCloseEvent>
#include <QtGui/QDesktopServices>
#include <QtGui/QDesktopWidget>
#include <QtGui/QDockWidget>
#include <QtGui/QFileDialog>
#include <QtGui/QHeaderView>
#include <QtGui/QInputDialog>
#include <QtGui/QListWidget>
#include <QtGui/QListWidgetItem>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QMessageBox>
#include <QtGui/QPainter>
#include <QtGui/QSplashScreen>
#include <QtGui/QStatusBar>
#include <QtGui/QTextEdit>
#include <QtGui/QToolBar>
#include <QtGui/QToolTip>
#include <QtGui/QToolButton>
#include <QtGui/QTreeWidget>
#include <QtGui/QTreeWidgetItem>
#include <QtGui/QWhatsThis>

#include <boost/math/special_functions/fpclassify.hpp>

#include <algorithm>
#include <utility>

using namespace std;

namespace OpenMS
{
  using namespace Internal;
	using namespace Math;

qreal TOPPViewBase::toppas_z_value_ = 42.0;
Int TOPPViewBase::toppas_node_offset_ = 0;

TOPPViewBase::TOPPViewBase(QWidget* parent):
        QMainWindow(parent),
        DefaultParamHandler("TOPPViewBase"),
        watcher_(0),
        watcher_msgbox_(false),
        toppas_clipboard_scene_(0)
        {
          setWindowTitle("TOPPView");
          setWindowIcon(QIcon(":/TOPPView.png"));

          //prevents errors caused by too small width,height values
          setMinimumSize(400,400);

          //enable drag-and-drop
          setAcceptDrops(true);

          //by default, linked zooming is turned off
          zoom_together_ = false;

          // get geometry of first screen
          QRect screen_geometry = QApplication::desktop()->screenGeometry();
          // center main window
          setGeometry(
            (int)(0.1 * screen_geometry.width()),
            (int)(0.1 * screen_geometry.height()),
            (int)(0.8 * screen_geometry.width()),
            (int)(0.8 * screen_geometry.height())
            );

          //set & create temporary path -- make sure its a new subdirectory, as TOPPASScene will delete it when its done
          toppas_tmp_path_ =  (File::getTempDirectory() + String(QDir::separator()) + File::getUniqueName()).toQString();
          QDir qd;
          qd.mkpath(toppas_tmp_path_);

          // create dummy widget (to be able to have a layout), Tab bar and workspace
          QWidget* dummy = new QWidget(this);
          setCentralWidget(dummy);
          QVBoxLayout* box_layout = new QVBoxLayout(dummy);

          // create empty tab bar and workspace which will hold the main visualization widgets (e.g. spectrawidgets...)
          tab_bar_ = new EnhancedTabBar(dummy);
          tab_bar_->setWhatsThis("Tab bar<BR><BR>Close tabs through the context menu or by double-clicking them.<BR>The tab bar accepts drag-and-drop from the layer bar.");
          tab_bar_->addTab("dummy",4710);
          tab_bar_->setMinimumSize(tab_bar_->sizeHint());
          tab_bar_->removeId(4710);

          connect(tab_bar_,SIGNAL(currentIdChanged(int)),this,SLOT(enhancedWorkspaceWindowChanged(int)));
          connect(tab_bar_,SIGNAL(aboutToCloseId(int)),this,SLOT(closeByTab(int)));

          //connect signals ans slots for drag-and-drop
          connect(tab_bar_,SIGNAL(dropOnWidget(const QMimeData*,QWidget*)),this,SLOT(copyLayer(const QMimeData*,QWidget*)));
          connect(tab_bar_,SIGNAL(dropOnTab(const QMimeData*,QWidget*,int)),this,SLOT(copyLayer(const QMimeData*,QWidget*,int)));
          box_layout->addWidget(tab_bar_);

          ws_= new EnhancedWorkspace(dummy);
          connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateToolBar()));
          connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateTabBar(QWidget*)));
          connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateLayerBar()));
          connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateViewBar()));
          connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateFilterBar()));
          connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateMenu()));
          connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateCurrentPath()));
          connect(ws_,SIGNAL(dropReceived(const QMimeData*,QWidget*,int)),this,SLOT(copyLayer(const QMimeData*,QWidget*,int)));

          box_layout->addWidget(ws_);

          //################## MENUS #################
          // File menu          
          QMenu* file = new QMenu("&File",this);
          menuBar()->addMenu(file);
          file->addAction("&Open file",this,SLOT(openFileDialog()), Qt::CTRL+Qt::Key_O);
          file->addAction("Open from &database",this,SLOT(openDatabaseDialog()), Qt::CTRL+Qt::Key_D);
          file->addAction("Open &example file",this,SLOT(openExampleDialog()));
          file->addAction("&Close",this,SLOT(closeFile()), Qt::CTRL+Qt::Key_W);
          file->addSeparator();

          // File menu: TOPPAS section
          file->addAction("&New TOPPAS pipeline",this,SLOT(newPipeline()), Qt::CTRL+Qt::Key_N);
          file->addAction("&Include TOPPAS pipeline",this,SLOT(includePipeline()), Qt::CTRL+Qt::Key_I);
          file->addAction("&Save TOPPAS pipeline",this,SLOT(savePipeline()), Qt::CTRL+Qt::Key_S);
          file->addAction("Save TOPPAS pipeline &As",this,SLOT(saveCurrentPipelineAs()), Qt::CTRL+Qt::SHIFT+Qt::Key_S);
          file->addAction("&Load TOPPAS resource file",this,SLOT(loadPipelineResourceFile()));
          file->addAction("Save TOPPAS &resource file",this,SLOT(savePipelineResourceFile()));
          file->addAction("Refresh TOPPAS &parameters",this,SLOT(refreshPipelineParameters()), Qt::CTRL+Qt::SHIFT+Qt::Key_P);
          file->addSeparator();

          //Meta data
          file->addAction("&Show meta data (file)",this,SLOT(metadataFileDialog()));
          file->addAction("&Show meta data (database)",this,SLOT(metadataDatabaseDialog()));
          file->addSeparator();

          //Recent files
          QMenu* recent_menu = new QMenu("&Recent files", this);
          recent_actions_.resize(20);
          for (Size i = 0; i<20; ++i)
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
          tools->addAction("&Go to",this,SLOT(showGoToDialog()), Qt::CTRL+Qt::Key_G);
          tools->addAction("&Edit meta data",this,SLOT(editMetadata()), Qt::CTRL+Qt::Key_M);
          tools->addAction("&Statistics",this,SLOT(layerStatistics()));
          tools->addSeparator();

          tools->addAction("Apply TOPP tool (whole layer)", this, SLOT(showTOPPDialog()), Qt::CTRL+Qt::Key_T)->setData(false);
          tools->addAction("Apply TOPP tool (visible layer data)", this, SLOT(showTOPPDialog()), Qt::CTRL+Qt::SHIFT+Qt::Key_T)->setData(true);
          tools->addAction("Rerun TOPP tool", this, SLOT(rerunTOPPTool()),Qt::Key_F4);
          tools->addSeparator();
          tools->addAction("&Annotate with identification", this, SLOT(annotateWithID()), Qt::CTRL+Qt::Key_I);
          tools->addAction("Align spectra", this, SLOT(showSpectrumAlignmentDialog()));
          tools->addAction("Generate theoretical spectrum", this, SLOT(showSpectrumGenerationDialog()));

          // TOPPAS menu
          // pipeline menu (TOPPAS pipeline)
          QMenu* toppas_pipeline = new QMenu("TOPPAS &Pipeline", this);
          toppas_pipeline->addAction("&Run (F5)",this,SLOT(runPipeline()));
          toppas_pipeline->addAction("&Abort",this,SLOT(abortPipeline()));
          menuBar()->addMenu(toppas_pipeline);

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
          QMenu* windows = new QMenu("&Windows", this);
          menuBar()->addMenu(windows);
          windows->addAction("&Cascade",this->ws_,SLOT(cascade()));
          windows->addAction("&Tile automatic",this->ws_,SLOT(tile()));
          windows->addAction(QIcon(":/tile_horizontal.png"),"Tile &vertical",this,SLOT(tileHorizontal()));
          windows->addAction(QIcon(":/tile_vertical.png"),"Tile &horizontal",this,SLOT(tileVertical()));
          linkZoom_action_ = windows->addAction("Link &Zoom",this,SLOT(linkZoom()));
          windows->addSeparator();

          //Help menu
          QMenu* help = new QMenu("&Help", this);
          menuBar()->addMenu(help);
          help->addAction(QWhatsThis::createAction(help));
          help->addSeparator();
          QAction* action = help->addAction("OpenMS website",this,SLOT(showURL()));
          action->setData("http://www.OpenMS.de");
          action = help->addAction("Tutorials and documentation",this,SLOT(showURL()), Qt::Key_F1);
          action->setData(String(File::getOpenMSDataPath() + "/../../doc/html/index.html").toQString());

          help->addSeparator();
          help->addAction("&About",this,SLOT(showAboutDialog()));

          //create status bar
          message_label_ = new QLabel(statusBar());
          statusBar()->addWidget(message_label_,1);

          rt_label_ = new QLabel("RT: 12345678", statusBar());
          rt_label_->setMinimumSize(rt_label_->sizeHint());
          rt_label_->setText("");
          statusBar()->addPermanentWidget(rt_label_,0);
          mz_label_ = new QLabel("m/z: 123456780912", statusBar());
          mz_label_->setMinimumSize(mz_label_->sizeHint());
          mz_label_->setText("");
          statusBar()->addPermanentWidget(mz_label_,0);

          //################## TOOLBARS #################
          //create toolbars and connect signals
          QToolButton* b;

          //--Basic tool bar for all views--
          tool_bar_ = addToolBar("Basic tool bar");

          //intensity modes
          intensity_button_group_ = new QButtonGroup(tool_bar_);
          intensity_button_group_->setExclusive(true);

          b = new QToolButton(tool_bar_);
          b->setIcon(QIcon(":/lin.png"));
          b->setToolTip("Intensity: Normal");
          b->setShortcut(Qt::Key_N);
          b->setCheckable(true);
          b->setWhatsThis("Intensity: Normal<BR><BR>Intensity is displayed unmodified.<BR>(Hotkey: N)");
          intensity_button_group_->addButton(b,SpectrumCanvas::IM_NONE);
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
          intensity_button_group_->addButton(b,SpectrumCanvas::IM_PERCENTAGE);
          tool_bar_->addWidget(b);

          b = new QToolButton(tool_bar_);
          b->setIcon(QIcon(":/snap.png"));
          b->setToolTip("Intensity: Snap to maximum displayed intensity");
          b->setShortcut(Qt::Key_S);
          b->setCheckable(true);
          b->setWhatsThis("Intensity: Snap to maximum displayed intensity<BR><BR> In this mode the"
                          " color gradient is adapted to the maximum currently displayed intensity."
                          "<BR>(Hotkey: S)");
          intensity_button_group_->addButton(b,SpectrumCanvas::IM_SNAP);
          tool_bar_->addWidget(b);

          b = new QToolButton(tool_bar_);
          b->setIcon(QIcon(":/log.png"));
          b->setToolTip("Intensity: Use log scaling for colors");
          b->setCheckable(true);
          b->setWhatsThis("Intensity: Logarithmic scaling of intensities for color calculation");
          intensity_button_group_->addButton(b,SpectrumCanvas::IM_LOG);
          tool_bar_->addWidget(b);

          connect(intensity_button_group_,SIGNAL(buttonClicked(int)),this,SLOT(setIntensityMode(int)));
          tool_bar_->addSeparator();

          //common buttons
          QAction* reset_zoom_button = tool_bar_->addAction(QIcon(":/reset_zoom.png"), "Reset Zoom", this, SLOT(resetZoom()));
          reset_zoom_button->setWhatsThis("Reset zoom: Zooms out as far as possible and resets the zoom history.<BR>(Hotkey: Backspace)");

          tool_bar_->show();

          //--1D toolbar--
          tool_bar_1d_ = addToolBar("1D tool bar");

          //draw modes 1D
          draw_group_1d_ = new QButtonGroup(tool_bar_1d_);
          draw_group_1d_->setExclusive(true);

          b = new QToolButton(tool_bar_1d_);
          b->setIcon(QIcon(":/peaks.png"));
          b->setToolTip("Peak mode");
          b->setShortcut(Qt::Key_I);
          b->setCheckable(true);
          b->setWhatsThis("1D Draw mode: Peaks<BR><BR>Peaks are diplayed as sticks.");
          draw_group_1d_->addButton(b,Spectrum1DCanvas::DM_PEAKS);
          tool_bar_1d_->addWidget(b);

          b = new QToolButton(tool_bar_1d_);
          b->setIcon(QIcon(":/lines.png"));
          b->setToolTip("Raw data mode");
          b->setShortcut(Qt::Key_R);
          b->setCheckable(true);
          b->setWhatsThis("1D Draw mode: Raw data<BR><BR>Peaks are diplayed as a continous line.");
          draw_group_1d_->addButton(b,Spectrum1DCanvas::DM_CONNECTEDLINES);
          tool_bar_1d_->addWidget(b);

          connect(draw_group_1d_,SIGNAL(buttonClicked(int)),this,SLOT(setDrawMode1D(int)));
          tool_bar_->addSeparator();

          //--2D peak toolbar--
          tool_bar_2d_peak_ = addToolBar("2D peak tool bar");

          dm_precursors_2d_ = tool_bar_2d_peak_->addAction(QIcon(":/precursors.png"),"Show fragment scan precursors");
          dm_precursors_2d_->setCheckable(true);
          dm_precursors_2d_->setWhatsThis("2D peak draw mode: Precursors<BR><BR>fragment scan precursor peaks are marked.<BR>(Hotkey: 1)");
          dm_precursors_2d_->setShortcut(Qt::Key_1);

          connect(dm_precursors_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

          projections_2d_ = tool_bar_2d_peak_->addAction(QIcon(":/projections.png"), "Show Projections" ,this, SLOT(toggleProjections()));
          projections_2d_->setWhatsThis("Projections: Shows projections of peak data along RT and MZ axis.<BR>(Hotkey: 2)");
          projections_2d_->setShortcut(Qt::Key_2);

          //--2D feature toolbar--
          tool_bar_2d_feat_ = addToolBar("2D feature tool bar");

          dm_hull_2d_ = tool_bar_2d_feat_->addAction(QIcon(":/convexhull.png"),"Show feature convex hull");
          dm_hull_2d_->setCheckable(true);
          dm_hull_2d_->setWhatsThis("2D feature draw mode: Convex hull<BR><BR>The convex hull of the feature is displayed.<BR>(Hotkey: 5)");
          dm_hull_2d_->setShortcut(Qt::Key_5);
          connect(dm_hull_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

          dm_hulls_2d_ = tool_bar_2d_feat_->addAction(QIcon(":/convexhulls.png"),"Show feature convex hulls");
          dm_hulls_2d_->setCheckable(true);
          dm_hulls_2d_->setWhatsThis("2D feature draw mode: Convex hulls<BR><BR>The convex hulls of the feature are displayed: One for each mass trace.<BR>(Hotkey: 6)");
          dm_hulls_2d_->setShortcut(Qt::Key_6);
          connect(dm_hulls_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

          // feature labels:
          dm_label_2d_ = new QToolButton(tool_bar_2d_feat_);
          dm_label_2d_->setPopupMode(QToolButton::MenuButtonPopup);
          QAction* action2 = new QAction(QIcon(":/labels.png"), "Show feature label", dm_label_2d_);
          action2->setCheckable(true);
          action2->setWhatsThis("2D feature draw mode: Labels<BR><BR>Display different kinds of annotation next to features.<BR>(Hotkey: 7)");
          action2->setShortcut(Qt::Key_7);
          dm_label_2d_->setDefaultAction(action2);
          tool_bar_2d_feat_->addWidget(dm_label_2d_);
          connect(dm_label_2d_, SIGNAL(triggered(QAction*)), this, SLOT(changeLabel(QAction*)));
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
          connect(dm_unassigned_2d_, SIGNAL(triggered(QAction*)), this, SLOT(changeUnassigned(QAction*)));
          //button menu
          group_unassigned_2d_ = new QActionGroup(dm_unassigned_2d_);
          menu = new QMenu(dm_unassigned_2d_);
          StringList options = StringList::create(
              "Don't show,Show by precursor m/z,Show by peptide mass");
          for (StringList::iterator opt_it = options.begin(); opt_it != options.end();
          ++opt_it)
          {
            QAction* temp = group_unassigned_2d_->addAction(opt_it->toQString());
            temp->setCheckable(true);
            if (opt_it == options.begin()) temp->setChecked(true);
            menu->addAction(temp);
          }
          dm_unassigned_2d_->setMenu(menu);

          //--2D consensus toolbar--
          tool_bar_2d_cons_ = addToolBar("2D peak tool bar");

          dm_elements_2d_ = tool_bar_2d_cons_->addAction(QIcon(":/elements.png"),"Show consensus feature element positions");
          dm_elements_2d_->setCheckable(true);
          dm_elements_2d_->setWhatsThis("2D consensus feature draw mode: Elements<BR><BR>The individual elements that make up the  consensus feature are drawn.<BR>(Hotkey: 9)");
          dm_elements_2d_->setShortcut(Qt::Key_9);
          connect(dm_elements_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

          //--2D identifications toolbar--
          tool_bar_2d_ident_ = addToolBar("2D identifications tool bar");

          dm_ident_2d_ = tool_bar_2d_ident_->addAction(QIcon(":/peptidemz.png"), "Use theoretical peptide mass for m/z positions (default: precursor mass)");
          dm_ident_2d_->setCheckable(true);
          dm_ident_2d_->setWhatsThis("2D peptide identification draw mode: m/z source<BR><BR>Toggle between precursor mass (default) and theoretical peptide mass as source for the m/z positions of peptide identifications.<BR>(Hotkey: 5)");
          dm_ident_2d_->setShortcut(Qt::Key_5);
          connect(dm_ident_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

          //################## Dock widgets #################
          // layer dock widget
          layer_dock_widget_ = new QDockWidget("Layers", this);
          addDockWidget(Qt::RightDockWidgetArea, layer_dock_widget_);
          layer_manager_ = new QListWidget(layer_dock_widget_);
          layer_manager_->setWhatsThis("Layer bar<BR><BR>Here the available layers are shown. Left-click on a layer to select it.<BR>Layers can be shown and hidden using the checkboxes in front of the name.<BR> Renaming and removing a layer is possible through the context menu.<BR>Dragging a layer to the tab bar copies the layer.<BR>Double-clicking a layer open its preferences.<BR>You can use the 'PageUp' and 'PageDown' buttons to change the selected layer.");

          layer_dock_widget_->setWidget(layer_manager_);
          layer_manager_->setContextMenuPolicy(Qt::CustomContextMenu);
          layer_manager_->setDragEnabled(true);
          connect(layer_manager_,SIGNAL(currentRowChanged(int)),this,SLOT(layerSelectionChange(int)));
          connect(layer_manager_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(layerContextMenu(const QPoint&)));
          connect(layer_manager_,SIGNAL(itemChanged(QListWidgetItem*)),this,SLOT(layerVisibilityChange(QListWidgetItem*)));
          connect(layer_manager_,SIGNAL(itemDoubleClicked(QListWidgetItem*)),this,SLOT(layerEdit(QListWidgetItem*)));

          windows->addAction(layer_dock_widget_->toggleViewAction());

          // Views dock widget
          views_dockwidget_ = new QDockWidget("Views", this);
          addDockWidget(Qt::RightDockWidgetArea, views_dockwidget_);
          views_tabwidget_ = new QTabWidget(views_dockwidget_);
          views_dockwidget_->setWidget(views_tabwidget_);

          spectra_view_widget_ = new SpectraViewWidget();
          connect(spectra_view_widget_, SIGNAL(showSpectrumMetaData(int)), this, SLOT(showSpectrumMetaData(int)));
          connect(spectra_view_widget_, SIGNAL(showSpectrumAs1D(int)), this, SLOT(showSpectrumAs1D(int)));
          connect(spectra_view_widget_, SIGNAL(showSpectrumAs1D(std::vector<int, std::allocator<int> >)), this, SLOT(showSpectrumAs1D(std::vector<int, std::allocator<int> >)));
          connect(spectra_view_widget_, SIGNAL(spectrumSelected(int)), this, SLOT(activate1DSpectrum(int)));
          connect(spectra_view_widget_, SIGNAL(spectrumSelected(std::vector<int, std::allocator<int> >)), this, SLOT(activate1DSpectrum(std::vector<int, std::allocator<int> >)));
          connect(spectra_view_widget_, SIGNAL(spectrumDoubleClicked(int)), this, SLOT(showSpectrumAs1D(int)));
          //connect(spectra_view_widget_, SIGNAL(spectrumDoubleClicked(int)), this, SLOT(showSpectrumAs1D(std::vector<int, std::allocator<int> >)));

          spectraview_behavior_ = new TOPPViewSpectraViewBehavior(this);
          view_behavior_ = spectraview_behavior_;

          spectra_identification_view_widget_ = new SpectraIdentificationViewWidget(Param());
          connect(spectra_identification_view_widget_, SIGNAL(spectrumDeselected(int)), this, SLOT(deactivate1DSpectrum(int)));
          connect(spectra_identification_view_widget_, SIGNAL(showSpectrumAs1D(int)), this, SLOT(showSpectrumAs1D(int)));
          // connect(spectra_identification_view_widget_, SIGNAL(showSpectrumAs1D(std::vector<int, std::allocator<int> >)), this, SLOT(showSpectrumAs1D(std::vector<int, std::allocator<int> >)));
          connect(spectra_identification_view_widget_, SIGNAL(spectrumSelected(int)), this, SLOT(activate1DSpectrum(int)));
          //connect(spectra_identification_view_widget_, SIGNAL(spectrumSelected(std::vector<int, std::allocator<int> >)), this, SLOT(activate1DSpectrum(std::vector<int, std::allocator<int> >)));
          identificationview_behavior_ = new TOPPViewIdentificationViewBehavior(this);
          connect(spectra_identification_view_widget_, SIGNAL(requestVisibleArea1D(DoubleReal, DoubleReal)), identificationview_behavior_, SLOT(setVisibleArea1D(DoubleReal, DoubleReal)));

          // topp tool dock widget          
          TOPPASTreeView* tools_tree_view = TOPPASBase::createTOPPToolsTreeWidget();
          connect (tools_tree_view, SIGNAL(itemDoubleClicked(QTreeWidgetItem*,int)), this, SLOT(insertNewVertexInCenter_(QTreeWidgetItem*)));

          views_tabwidget_->addTab(spectra_view_widget_, "Scan view");
          views_tabwidget_->addTab(spectra_identification_view_widget_, "Identification view");
          views_tabwidget_->setTabEnabled(1, false);
          views_tabwidget_->addTab(tools_tree_view, "TOPPAS view");
          views_tabwidget_->setTabEnabled(2, false);

          // switch between different view tabs
          connect(views_tabwidget_, SIGNAL(currentChanged(int)), this, SLOT(viewChanged(int)));

          // add hide/show option to dock widget
          windows->addAction(views_dockwidget_->toggleViewAction());

          // filter dock widget
          filter_dock_widget_ = new QDockWidget("Data filters", this);
          addDockWidget(Qt::RightDockWidgetArea, filter_dock_widget_);
          QWidget* tmp_widget = new QWidget(); //dummy widget as QDockWidget takes only one widget
          filter_dock_widget_->setWidget(tmp_widget);

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

          windows->addAction(filter_dock_widget_->toggleViewAction());

          //log window
          QDockWidget* log_bar = new QDockWidget("Log", this);
          addDockWidget(Qt::BottomDockWidgetArea, log_bar);
          log_ = new QTextEdit(log_bar);
          log_->setReadOnly(true);
          log_->setContextMenuPolicy(Qt::CustomContextMenu);
          connect(log_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(logContextMenu(const QPoint&)));
          log_bar->setWidget(log_);
          log_bar->hide();
          windows->addAction(log_bar->toggleViewAction());

          //################## DEFAULTS #################
          initializeDefaultParameters_();

          // store defaults in param_
          defaultsToParam_();

          //load param file
          loadPreferences();

          //set current path
          current_path_ = param_.getValue("preferences:default_path");

          //update the menu
          updateMenu();

          topp_.process = 0;

          //######################### File System Watcher ###########################################
          watcher_ = new FileWatcher(this);
          connect(watcher_,SIGNAL(fileChanged(const String&)),this, SLOT(fileChanged_(const String&)));
        }

  void TOPPViewBase::initializeDefaultParameters_()
  {
    //general
    defaults_.setValue("preferences:default_map_view", "2d", "Default visualization mode for maps.");
    defaults_.setValidStrings("preferences:default_map_view",StringList::create("2d,3d"));
    defaults_.setValue("preferences:default_path", ".", "Default path for loading and storing files.");
    defaults_.setValue("preferences:default_path_current", "true", "If the current path is preferred over the default path.");
    defaults_.setValidStrings("preferences:default_path_current",StringList::create("true,false"));
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
    defaults_.setValue("preferences:topp_cleanup", "true", "If the temporary files for calling of TOPP tools should be removed after the call.");
    defaults_.setValidStrings("preferences:topp_cleanup",StringList::create("true,false"));
    //db
    defaults_.setValue("preferences:db:host", "localhost", "Database server host name.");
    defaults_.setValue("preferences:db:login", "NoName", "Database login.");
    defaults_.setValue("preferences:db:name", "OpenMS", "Database name.");
    defaults_.setValue("preferences:db:port", 3306, "Database server port.");
    defaults_.setSectionDescription("preferences:db","Database settings.");
    // 1d view
    Spectrum1DCanvas* def1 = new Spectrum1DCanvas(Param(),0);
    defaults_.insert("preferences:1d:",def1->getDefaults());
    delete def1;
    defaults_.setSectionDescription("preferences:1d","Settings for single spectrum view.");
    // 2d view
    Spectrum2DCanvas* def2 = new Spectrum2DCanvas(Param(),0);
    defaults_.insert("preferences:2d:",def2->getDefaults());
    defaults_.setSectionDescription("preferences:2d","Settings for 2D map view.");
    delete def2;
    // 3d view
    Spectrum3DCanvas* def3 = new Spectrum3DCanvas(Param(),0);
    defaults_.insert("preferences:3d:",def3->getDefaults());
    delete def3;
    defaults_.setSectionDescription("preferences:3d","Settings for 3D map view.");
    // identification view
    SpectraIdentificationViewWidget* def4 = new SpectraIdentificationViewWidget(Param(),0);
    defaults_.insert("preferences:idview:",def4->getDefaults());
    delete def4;
    defaults_.setSectionDescription("preferences:idview","Settings for identification view.");
    defaults_.setValue("preferences:version","none","OpenMS version, used to check if the TOPPView.ini is up-to-date");
    subsections_.push_back("preferences:RecentFiles");
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

    // create managed pointer to experiment data
    ExperimentType* exp = new ExperimentType();
    ExperimentSharedPtrType exp_sptr(exp);

    FeatureMapType* dummy_map = new FeatureMapType();
    FeatureMapSharedPtrType dummy_map_sptr(dummy_map);

    ConsensusMapType* dummy_map2 = new ConsensusMapType();
    ConsensusMapSharedPtrType dummy_map2_sptr(dummy_map2);

    vector<PeptideIdentification> dummy_peptides;
    try
    {
      db.loadExperiment(db_id, *exp);
    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::critical(this,"Error",(String("Error while reading data: ")+e.what()).c_str());
      setCursor(Qt::ArrowCursor);
      return;
    }
    exp_sptr->sortSpectra(true);
    exp_sptr->updateRanges(1);

    //determine if the data is 1D or 2D
    QSqlQuery result = con.executeQuery(String("SELECT count(id) from DATA_Spectrum where fid_MSExperiment='")+db_id+"' and MSLevel='1'");
    LayerData::DataType data_type = ((result.value(0).toInt()>1) ?
                                     LayerData::DT_PEAK :
                                     LayerData::DT_CHROMATOGRAM);

    //add data
    if (caption=="") caption = String("DB entry ")+db_id;
    addData(dummy_map_sptr, dummy_map2_sptr, dummy_peptides, exp_sptr, data_type, false, show_options, true, "", caption, window_id);

    //Reset cursor
    setCursor(Qt::ArrowCursor);
  }

  // static
  bool TOPPViewBase::containsMS1Scans(const ExperimentType& exp)
  {
    //test if no scans with MS-level 1 exist => prevent deadlock
    bool ms1_present = false;
    for (Size i = 0; i < exp.size(); ++i)
    {
        if (exp[i].getMSLevel() == 1)
        {
            ms1_present = true;
            break;
        }
    }
    return ms1_present;
  }

  // static
  float TOPPViewBase::estimateNoiseFromRandomMS1Scans(const ExperimentType& exp, UInt n_scans)
  {
    if (!TOPPViewBase::containsMS1Scans(exp))
    {
      return 0.0;
    }

    float noise = 0.0;
    UInt count = 0;
    srand(time(0));
    while (count < n_scans)
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
        noise += tmp[(UInt)ceil((float)(tmp.size()-1)/1.25f)];
        ++count;
      }
    }
    return noise / (DoubleReal)n_scans;
  }

  // static
  UInt TOPPViewBase::countMS1Zeros(const ExperimentType& exp)
  {
    if (!TOPPViewBase::containsMS1Scans(exp))
    {
      return 0;
    }

    UInt zeros = 0;
    for (Size i = 0; i != exp.size(); ++i)
    {
      if (exp[i].getMSLevel() != 1) // skip non MS1-level scans
      {
        continue;
      }
      for(Size j = 0; j != exp[i].size(); ++j)
      {
        DoubleReal intensity = exp[i][j].getIntensity();
        if (intensity == 0.0)
        {
          zeros++;
        }
      }
    }
    return zeros;
  }

  // static
  bool TOPPViewBase::hasPeptideIdentifications(const ExperimentType& map)
  {
    for(Size i = 0; i!= map.size(); ++i)
    {
      if ( !map[i].getPeptideIdentifications().empty())
      {
        return true;
      }
    }
    return false;
  }

  void TOPPViewBase::preferencesDialog()
  {
		Internal::TOPPViewPrefDialog dlg(this);

    // --------------------------------------------------------------------
    // Get pointers to the widget in the preferences dialog

    // default tab
		QLineEdit* default_path = dlg.findChild<QLineEdit*>("default_path");
		QCheckBox* default_path_current = dlg.findChild<QCheckBox*>("default_path_current");
		QLineEdit* temp_path = dlg.findChild<QLineEdit*>("temp_path");
		QSpinBox* recent_files = dlg.findChild<QSpinBox*>("recent_files");
		QComboBox* map_default = dlg.findChild<QComboBox*>("map_default");
		QComboBox* map_cutoff = dlg.findChild<QComboBox*>("map_cutoff");
		QComboBox* on_file_change = dlg.findChild<QComboBox*>("on_file_change");

    // db tab
		QLineEdit* db_host = dlg.findChild<QLineEdit*>("db_host");
		QSpinBox* db_port = dlg.findChild<QSpinBox*>("db_port");
		QLineEdit* db_name = dlg.findChild<QLineEdit*>("db_name");
		QLineEdit* db_login = dlg.findChild<QLineEdit*>("db_login");

    // 1D view tab
		ColorSelector* color_1D = dlg.findChild<ColorSelector*>("color_1D");
		ColorSelector* selected_1D = dlg.findChild<ColorSelector*>("selected_1D");
		ColorSelector* icon_1D = dlg.findChild<ColorSelector*>("icon_1D");

    // 2D view tab
		MultiGradientSelector* peak_2D = dlg.findChild<MultiGradientSelector*>("peak_2D");
		QComboBox* mapping_2D = dlg.findChild<QComboBox*>("mapping_2D");
		QComboBox* feature_icon_2D = dlg.findChild<QComboBox*>("feature_icon_2D");
		QSpinBox* feature_icon_size_2D = dlg.findChild<QSpinBox*>("feature_icon_size_2D");

    // 3D view tab
		MultiGradientSelector* peak_3D = dlg.findChild<MultiGradientSelector*>("peak_3D");
		QComboBox* shade_3D = dlg.findChild<QComboBox*>("shade_3D");
		QSpinBox* line_width_3D  = dlg.findChild<QSpinBox*>("line_width_3D");

    // identification view tab
    QListWidget* id_view_ions = dlg.findChild<QListWidget*>("ions_list_widget");
    QDoubleSpinBox* a_intensity = dlg.findChild<QDoubleSpinBox*>("a_intensity");
    QDoubleSpinBox* b_intensity = dlg.findChild<QDoubleSpinBox*>("b_intensity");
    QDoubleSpinBox* c_intensity = dlg.findChild<QDoubleSpinBox*>("c_intensity");
    QDoubleSpinBox* x_intensity = dlg.findChild<QDoubleSpinBox*>("x_intensity");
    QDoubleSpinBox* y_intensity = dlg.findChild<QDoubleSpinBox*>("y_intensity");
    QDoubleSpinBox* z_intensity = dlg.findChild<QDoubleSpinBox*>("z_intensity");

    QDoubleSpinBox* tolerance = dlg.findChild<QDoubleSpinBox*>("tolerance");
    QCheckBox* is_relative_tolerance = dlg.findChild<QCheckBox*>("unit_is_ppm");

    QDoubleSpinBox* relative_loss_intensity = dlg.findChild<QDoubleSpinBox*>("relative_loss_intensity");
    QSpinBox* max_isotopes = dlg.findChild<QSpinBox*>("max_isotopes");
    QSpinBox* charge = dlg.findChild<QSpinBox*>("charge");

    QList<QListWidgetItem*> a_ions = id_view_ions->findItems("A-ions", Qt::MatchFixedString);
    QList<QListWidgetItem*> b_ions = id_view_ions->findItems("B-ions", Qt::MatchFixedString);
    QList<QListWidgetItem*> c_ions = id_view_ions->findItems("C-ions", Qt::MatchFixedString);
    QList<QListWidgetItem*> x_ions = id_view_ions->findItems("X-ions", Qt::MatchFixedString);
    QList<QListWidgetItem*> y_ions = id_view_ions->findItems("Y-ions", Qt::MatchFixedString);
    QList<QListWidgetItem*> z_ions = id_view_ions->findItems("Z-ions", Qt::MatchFixedString);
    QList<QListWidgetItem*> pc_ions = id_view_ions->findItems("Precursor", Qt::MatchFixedString);
    QList<QListWidgetItem*> nl_ions = id_view_ions->findItems("Neutral losses", Qt::MatchFixedString);
    QList<QListWidgetItem*> ic_ions = id_view_ions->findItems("Isotope clusters", Qt::MatchFixedString);
    QList<QListWidgetItem*> ai_ions = id_view_ions->findItems("Abundant immonium-ions", Qt::MatchFixedString);

    // --------------------------------------------------------------------
    // Set dialog entries from current parameter object (default values)

    // default
		default_path->setText(param_.getValue("preferences:default_path").toQString());
		if ((String)param_.getValue("preferences:default_path_current")=="true")
		{
			default_path_current->setChecked(true);
		}
		else
		{
			default_path_current->setChecked(false);
		}
		temp_path->setText(param_.getValue("preferences:tmp_file_path").toQString());
		recent_files->setValue((Int)param_.getValue("preferences:number_of_recent_files"));
		map_default->setCurrentIndex(map_default->findText(param_.getValue("preferences:default_map_view").toQString()));
		map_cutoff->setCurrentIndex(map_cutoff->findText(param_.getValue("preferences:intensity_cutoff").toQString()));
		on_file_change->setCurrentIndex(on_file_change->findText(param_.getValue("preferences:on_file_change").toQString()));

    // db
		db_host->setText(param_.getValue("preferences:db:host").toQString());
		db_port->setValue((Int)param_.getValue("preferences:db:port"));
		db_name->setText(param_.getValue("preferences:db:name").toQString());
		db_login->setText(param_.getValue("preferences:db:login").toQString());

    // 1D view
		color_1D->setColor(QColor(param_.getValue("preferences:1d:peak_color").toQString()));
		selected_1D->setColor(QColor(param_.getValue("preferences:1d:highlighted_peak_color").toQString()));
		icon_1D->setColor(QColor(param_.getValue("preferences:1d:icon_color").toQString()));

    // 2D view
		peak_2D->gradient().fromString(param_.getValue("preferences:2d:dot:gradient"));
		mapping_2D->setCurrentIndex(mapping_2D->findText(param_.getValue("preferences:2d:mapping_of_mz_to").toQString()));
		feature_icon_2D->setCurrentIndex(feature_icon_2D->findText(param_.getValue("preferences:2d:dot:feature_icon").toQString()));
		feature_icon_size_2D->setValue((Int)param_.getValue("preferences:2d:dot:feature_icon_size"));

    // 3D view
		peak_3D->gradient().fromString(param_.getValue("preferences:3d:dot:gradient"));
		shade_3D->setCurrentIndex((Int)param_.getValue("preferences:3d:dot:shade_mode"));
		line_width_3D->setValue((Int)param_.getValue("preferences:3d:dot:line_width"));

    // id view
    a_intensity->setValue((DoubleReal)param_.getValue("preferences:idview:a_intensity"));
    b_intensity->setValue((DoubleReal)param_.getValue("preferences:idview:b_intensity"));
    c_intensity->setValue((DoubleReal)param_.getValue("preferences:idview:c_intensity"));
    x_intensity->setValue((DoubleReal)param_.getValue("preferences:idview:x_intensity"));
    y_intensity->setValue((DoubleReal)param_.getValue("preferences:idview:y_intensity"));
    z_intensity->setValue((DoubleReal)param_.getValue("preferences:idview:z_intensity"));
    tolerance->setValue((DoubleReal)param_.getValue("preferences:idview:tolerance"));

    relative_loss_intensity->setValue((DoubleReal)param_.getValue("preferences:idview:relative_loss_intensity"));
    max_isotopes->setValue((Int)param_.getValue("preferences:idview:max_isotope"));
    charge->setValue((Int)param_.getValue("preferences:idview:charge"));

    if(a_ions.empty())
    {
      showLogMessage_(LS_ERROR,"", "String 'A-ions' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:show_a_ions").toBool() == true ? Qt::Checked : Qt::Unchecked;
      a_ions[0]->setCheckState(state);
    }

    if(b_ions.empty())
    {
      showLogMessage_(LS_ERROR,"", "String 'B-ions' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:show_b_ions").toBool() == true ? Qt::Checked : Qt::Unchecked;
      b_ions[0]->setCheckState(state);
    }

    if(c_ions.empty())
    {
      showLogMessage_(LS_ERROR,"", "String 'C-ions' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:show_c_ions").toBool() == true ? Qt::Checked : Qt::Unchecked;
      c_ions[0]->setCheckState(state);
    }

    if(x_ions.empty())
    {
      showLogMessage_(LS_ERROR,"", "String 'X-ions' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:show_x_ions").toBool() == true ? Qt::Checked : Qt::Unchecked;
      x_ions[0]->setCheckState(state);
    }

    if(y_ions.empty())
    {
      showLogMessage_(LS_ERROR,"", "String 'Y-ions' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:show_y_ions").toBool() == true ? Qt::Checked : Qt::Unchecked;
      y_ions[0]->setCheckState(state);
    }

    if(z_ions.empty())
    {
      showLogMessage_(LS_ERROR,"", "String 'Z-ions' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:show_z_ions").toBool() == true ? Qt::Checked : Qt::Unchecked;
      z_ions[0]->setCheckState(state);
    }

    if(pc_ions.empty())
    {
      showLogMessage_(LS_ERROR,"", "String 'Precursor' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:show_precursor").toBool() == true ? Qt::Checked : Qt::Unchecked;
      pc_ions[0]->setCheckState(state);
    }

    if(nl_ions.empty())
    {
      showLogMessage_(LS_ERROR,"", "String 'Neutral losses' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:add_losses").toBool() == true ? Qt::Checked : Qt::Unchecked;
      nl_ions[0]->setCheckState(state);
    }

    if(ic_ions.empty())
    {
      showLogMessage_(LS_ERROR,"", "String 'Isotope clusters' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:add_isotopes").toBool() == true ? Qt::Checked : Qt::Unchecked;
      ic_ions[0]->setCheckState(state);
    }

    if(ai_ions.empty())
    {
      showLogMessage_(LS_ERROR,"", "String 'Abundant immonium-ions' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:add_abundant_immonium_ions").toBool() == true ? Qt::Checked : Qt::Unchecked;
      ai_ions[0]->setCheckState(state);
    }

    if(is_relative_tolerance == 0)
    {
      showLogMessage_(LS_ERROR,"", "String 'unit is ppm' doesn't exist in identification dialog.");
    } else
    {
      Qt::CheckState state = param_.getValue("preferences:idview:is_relative_tolerance").toBool() == true ? Qt::Checked : Qt::Unchecked;
      is_relative_tolerance->setCheckState(state);
    }

    // --------------------------------------------------------------------
    // Execute dialog and update parameter object with user modified values
		if (dlg.exec())
		{
			param_.setValue("preferences:default_path", default_path->text());
			if (default_path_current->isChecked())
			{
				param_.setValue("preferences:default_path_current","true");
			}
			else
			{
				param_.setValue("preferences:default_path_current","false");
			}
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
			param_.setValue("preferences:2d:dot:feature_icon",feature_icon_2D->currentText());
			param_.setValue("preferences:2d:dot:feature_icon_size",feature_icon_size_2D->value());

			param_.setValue("preferences:3d:dot:gradient",peak_3D->gradient().toString());
			param_.setValue("preferences:3d:dot:shade_mode", shade_3D->currentIndex());
			param_.setValue("preferences:3d:dot:line_width",line_width_3D->value());

      // id view
      param_.setValue("preferences:idview:a_intensity", a_intensity->value(), "Default intensity of a-ions");
      param_.setValue("preferences:idview:b_intensity", b_intensity->value(), "Default intensity of b-ions");
      param_.setValue("preferences:idview:c_intensity", c_intensity->value(), "Default intensity of c-ions");
      param_.setValue("preferences:idview:x_intensity", x_intensity->value(), "Default intensity of x-ions");
      param_.setValue("preferences:idview:y_intensity", y_intensity->value(), "Default intensity of y-ions");
      param_.setValue("preferences:idview:z_intensity", z_intensity->value(), "Default intensity of z-ions");
      param_.setValue("preferences:idview:relative_loss_intensity", relative_loss_intensity->value(), "Relativ loss in percent");
      param_.setValue("preferences:idview:max_isotope", max_isotopes->value(), "Maximum number of isotopes");
      param_.setValue("preferences:idview:charge", charge->value(), "Charge state");
      param_.setValue("preferences:idview:tolerance", tolerance->value(), "Alignment tolerance");

      String checked;
      a_ions[0]->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:show_a_ions", checked, "Show a-ions");

      b_ions[0]->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:show_b_ions", checked, "Show b-ions");

      c_ions[0]->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:show_c_ions", checked, "Show c-ions");

      x_ions[0]->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:show_x_ions", checked, "Show x-ions");

      y_ions[0]->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:show_y_ions", checked, "Show y-ions");

      z_ions[0]->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:show_z_ions", checked, "Show z-ions");

      pc_ions[0]->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:show_precursor", checked, "Show precursor");

      nl_ions[0]->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:add_losses", checked, "Show neutral losses");

      ic_ions[0]->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:add_isotopes", checked, "Show isotopes");

      ai_ions[0]->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:add_abundant_immonium_ions", checked, "Show abundant immonium ions");

      is_relative_tolerance->checkState() == Qt::Checked ? checked = "true" : checked = "false";
      param_.setValue("preferences:idview:is_relative_tolerance", checked, "Use ppm instead of Da for the automatic alignment");
			savePreferences();
		}
  }

  std::set<String> TOPPViewBase::getFilenamesOfOpenFiles_()
  {
    set<String> filename_set;
    // iterate over all windows
    QWidgetList wl = ws_->windowList();
    for(int i=0; i!=ws_->windowList().count(); ++i)
    {
      QWidget* w = wl[i];
      // iterate over all widgets
      const SpectrumWidget* sw = qobject_cast<const SpectrumWidget*>(w);
      if (sw!=0)
      {
        Size lc = sw->canvas()->getLayerCount();
        // iterate over all layers
        for (Size j=0; j!= lc; ++j)
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

  	String abs_filename = File::absolutePath(filename);

    // check if the file exists
    if (!File::exists(abs_filename))
    {
      showLogMessage_(LS_ERROR, "Open file error",String("The file '") + abs_filename + "' does not exist!");
    	setCursor(Qt::ArrowCursor);
      return;
    }

    // determine file type
  	FileHandler fh;
		FileTypes::Type file_type = fh.getType(abs_filename);
    if (file_type == FileTypes::UNKNOWN)
		{
      showLogMessage_(LS_ERROR,"Open file error",String("Could not determine file type of '") + abs_filename + "'!");
    	setCursor(Qt::ArrowCursor);
      return;
		}

    // abort if file type unsupported
    if (file_type == FileTypes::INI)
		{
      showLogMessage_(LS_ERROR,"Open file error",String("The type '") + fh.typeToName(file_type) + "' is not supported!");
   		setCursor(Qt::ArrowCursor);
      return;
		}

		//try to load data and determine if it's 1D or 2D data

    // create shared pointer to main data types
    FeatureMapType* feature_map = new FeatureMapType();
    FeatureMapSharedPtrType feature_map_sptr(feature_map);

    ExperimentType* peak_map = new ExperimentType();
    ExperimentSharedPtrType peak_map_sptr(peak_map);

    ConsensusMapType* consensus_map = new ConsensusMapType();
    ConsensusMapSharedPtrType consensus_map_sptr(consensus_map);

		vector<PeptideIdentification> peptides;

		LayerData::DataType data_type;

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
				data_type = LayerData::DT_IDENT;
			}
      else
      {
        fh.loadExperiment(abs_filename, *peak_map, file_type, ProgressLogger::GUI);
        data_type = LayerData::DT_CHROMATOGRAM;
        if (TOPPViewBase::containsMS1Scans(*peak_map))
        {
          data_type = LayerData::DT_PEAK;
        }
      }
    }
    catch(Exception::BaseException& e)
    {
      showLogMessage_(LS_ERROR,"Error while loading file", e.what());
    	setCursor(Qt::ArrowCursor);
      return;
    }

    // sort for mz and update ranges of newly loaded data
    peak_map_sptr->sortSpectra(true);
    peak_map_sptr->updateRanges(1);

    // try to add the data
    if (caption == "")
		{
			caption = File::basename(abs_filename);
    }
    else
    {
    	abs_filename = "";
    }

    addData(feature_map_sptr, consensus_map_sptr, peptides, peak_map_sptr, data_type, false, show_options, true, abs_filename, caption, window_id, spectrum_id);

    // add to recent file
    if (add_to_recent)
    {
      addRecentFile_(filename);
    }

    // watch file contents for changes
    watcher_->addFile(abs_filename);

    // reset cursor
    setCursor(Qt::ArrowCursor);
  }

  void TOPPViewBase::addData(FeatureMapSharedPtrType feature_map, ConsensusMapSharedPtrType consensus_map, vector<PeptideIdentification>& peptides, ExperimentSharedPtrType peak_map, LayerData::DataType data_type, bool show_as_1d, bool show_options, bool as_new_window, const String& filename, const String& caption, UInt window_id, Size spectrum_id)
  {
    // initialize flags with defaults from the parameters
  	bool maps_as_2d = ((String)param_.getValue("preferences:default_map_view")=="2d");
  	bool maps_as_1d = false;
    bool use_intensity_cutoff = ((String)param_.getValue("preferences:intensity_cutoff")=="on");

    // feature, consensus feature and identifications can be merged
		bool mergeable = ((data_type == LayerData::DT_FEATURE) ||
											(data_type == LayerData::DT_CONSENSUS) ||
											(data_type == LayerData::DT_IDENT));

    // only one peak spectrum? disable 2D as default
    if (peak_map->size() == 1)
    {
      maps_as_2d = false;
    }

    // set the window where (new layer) data could be opened in
    // get EnhancedTabBarWidget with given id
    EnhancedTabBarWidgetInterface* tab_bar_target = window_(window_id);

    // cast to SpectrumWidget
    SpectrumWidget* target_window = dynamic_cast<SpectrumWidget*>(tab_bar_target);

    if (tab_bar_target == 0)
		{
      target_window = getActiveSpectrumWidget();
		}
    else
		{
			as_new_window = false;
		}

    //create dialog no matter if it is shown or not. It is used to determine the flags.
    TOPPViewOpenDialog dialog(caption, as_new_window, maps_as_2d, use_intensity_cutoff, this);

    //disable opening in new window when there is no active window or feature/ID data is to be opened, but the current window is a 3D window
    if (target_window == 0 || (mergeable && dynamic_cast<Spectrum3DWidget*>(target_window) != 0))
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
    if (mergeable && target_window != 0) //TODO merge
		{
      SpectrumCanvas* open_canvas = target_window->canvas();
			Map<Size,String> layers;
			for (Size i=0; i<open_canvas->getLayerCount(); ++i)
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
  	Int merge_layer = dialog.getMergeLayer();

		//determine the window to open the data in
		if (as_new_window) //new window
    {
      if (maps_as_1d) // 2d in 1d window
      {
        target_window = new Spectrum1DWidget(getSpectrumParameters(1), ws_);
      }
      else if (maps_as_2d || mergeable) //2d or features/IDs
      {
        target_window = new Spectrum2DWidget(getSpectrumParameters(2), ws_);
      }
      else // 3d
      {
        target_window = new Spectrum3DWidget(getSpectrumParameters(3), ws_);
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
        if (!target_window->canvas()->addLayer(consensus_map,filename)) return;
			}
			else if (data_type == LayerData::DT_IDENT)
			{
        if (!target_window->canvas()->addLayer(peptides, filename)) return;
			}
	    else //peaks
	    {
        if (!target_window->canvas()->addLayer(peak_map,filename)) return;

	      //calculate noise
        if (use_intensity_cutoff)
	      {
          DoubleReal cutoff = estimateNoiseFromRandomMS1Scans(*(target_window->canvas()->getCurrentLayer().getPeakData()));
					//create filter
					DataFilters::DataFilter filter;
					filter.field = DataFilters::INTENSITY;
					filter.op = DataFilters::GREATER_EQUAL;
					filter.value = cutoff;
					///add filter
					DataFilters filters;
					filters.add(filter);
          target_window->canvas()->setFilters(filters);
        } else  // no mower, hide zeros if wanted
        {
          Int n_zeros = TOPPViewBase::countMS1Zeros(*(target_window->canvas()->getCurrentLayer().getPeakData()));
          if (n_zeros > 0)
          {
            //create filter
            DataFilters::DataFilter filter;
            filter.field = DataFilters::INTENSITY;
            filter.op = DataFilters::GREATER_EQUAL;
            filter.value = 0.001;
            QMessageBox::question(this, "Note:", "Data contains zero values.\nA filter will be added to hide these values.\nYou can reenable data points with zero intensity by removing the filter.");
            ///add filter
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

			//set caption
      target_window->canvas()->setLayerName(target_window->canvas()->activeLayerIndex(), caption);
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
      showSpectrumWidgetInWindow(target_window,caption);
		}
    //updateDataBar();
		updateLayerBar();
    updateViewBar();
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
    // get/correct number of recent files
		UInt number_of_recent_files = UInt(param_.getValue("preferences:number_of_recent_files"));
		if (number_of_recent_files>20)
		{
			number_of_recent_files = 20;
			param_.setValue("preferences:number_of_recent_files",20);
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

  EnhancedTabBarWidgetInterface* TOPPViewBase::window_(int id) const
  {
    // return window with window_id == id
  	QList<QWidget*> windows = ws_->windowList();

    for(int i = 0; i < windows.size(); ++i)
    {
      EnhancedTabBarWidgetInterface* w = dynamic_cast<EnhancedTabBarWidgetInterface*>(windows.at(i));
      if (w->getWindowId() == id)
			{
        return w;
			}
		}
		return 0;
  }

  void TOPPViewBase::closeByTab(int id)
  {
    QWidget* w = dynamic_cast<QWidget*>(window_(id));
    if (w)
  	{
      w->close();
  		updateMenu();
  	}
  }

  void TOPPViewBase::enhancedWorkspaceWindowChanged(int id)
  {    
    QWidget* w = dynamic_cast<QWidget*>(window_(id));
    if (w)
  	{
      w->setFocus();
      TOPPASWidget* tw = dynamic_cast<TOPPASWidget*>(w);
      SpectrumWidget* sw = dynamic_cast<SpectrumWidget*>(w);
      if (tw)  // TOPPASWidget
      {
        setTOPPASTabEnabled(true);
        views_tabwidget_->setCurrentIndex(2);
        views_tabwidget_->setTabEnabled(0, false);  // switch scan view off
        views_tabwidget_->setTabEnabled(1, false);  // switch identification view off
      }
      else if (sw)  // SpectrumWidget
      {
        views_tabwidget_->setTabEnabled(0, true);

        // check if there is a layer before requesting data from it
        if(sw->canvas()->getLayerCount() > 0)
        {
          const ExperimentType& map = *sw->canvas()->getCurrentLayer().getPeakData();
          if(hasPeptideIdentifications(map))
          {
            views_tabwidget_->setTabEnabled(1, true);
            if (dynamic_cast<Spectrum2DWidget*>(w))
            {
              views_tabwidget_->setCurrentIndex(0);  // switch to scan tab for 2D widget
            } else if (dynamic_cast<Spectrum1DWidget*>(w))
            {
              views_tabwidget_->setCurrentIndex(1);  // switch to identification tab for 1D widget
            }
          }
          else
          {
            views_tabwidget_->setTabEnabled(1, false);
            views_tabwidget_->setCurrentIndex(0); // stay on scan view tab
          }
          setTOPPASTabEnabled(false);
        }
      }
  	}
  }

  void TOPPViewBase::closeFile()
  {
    ws_->activeWindow()->close();
    updateMenu();
  }

  void TOPPViewBase::editMetadata()
  {
    SpectrumCanvas* canvas = getActiveCanvas();

    // warn if hidden layer => wrong layer selected...
  	if (!canvas->getCurrentLayer().visible)
  	{
  		showLogMessage_(LS_NOTICE,"The current layer is not visible","Have you selected the right layer for this action?");
  	}

  	//show editable meta data dialog
  	canvas->showMetaData(true);
  }

  void TOPPViewBase::layerStatistics()
  {
    getActiveSpectrumWidget()->showStatistics();
  	updateFilterBar();
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
    QApplication::processEvents();
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
      mz_label_->setText((String("m/z: ")+String::number(mz,6).fillLeft(' ',8)).toQString());
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
      rt_label_->setText((String("RT: ")+String::number(rt,1).fillLeft(' ',8)).toQString());
    }
    statusBar()->update();
  }

  void TOPPViewBase::resetZoom()
  {
    SpectrumWidget* w = getActiveSpectrumWidget();
    if (w != 0)
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
      Spectrum2DWidget* w2d = dynamic_cast<Spectrum2DWidget*>(w);
      // 2D widget and intensity mode changed?
      if (w2d && w2d->canvas()->getIntensityMode() != index)
      {
        if (index == OpenMS::SpectrumCanvas::IM_LOG)
        {
          w2d->canvas()->getCurrentLayer().param.setValue("dot:gradient", MultiGradient::getDefaultGradientLogarithmicIntensityMode().toString());
          w2d->canvas()->recalculateCurrentLayerDotGradient();
        } else if (index != OpenMS::SpectrumCanvas::IM_LOG)
        {
          w2d->canvas()->getCurrentLayer().param.setValue("dot:gradient", MultiGradient::getDefaultGradientLinearIntensityMode().toString());
          w2d->canvas()->recalculateCurrentLayerDotGradient();
        }
      }
    	w->setIntensityMode((OpenMS::SpectrumCanvas::IntensityModes)index);
  	}
  }

  void TOPPViewBase::setDrawMode1D(int index)
  {
    Spectrum1DWidget* w = getActive1DWidget();
    if (w)
    {
      draw_group_1d_->button(Spectrum1DCanvas::DM_PEAKS)->setChecked(true);
    	w->canvas()->setDrawMode((OpenMS::Spectrum1DCanvas::DrawModes)index);
  	}
  }

	void TOPPViewBase::changeLabel(QAction* action)
	{
    bool set = false;

		//label type is selected
		for (Size i=0; i<LayerData::SIZE_OF_LABEL_TYPE; ++i)
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
		bool set = false;

		// mass reference is selected
		if (action->text().toStdString() == "Don't show")
		{
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::F_UNASSIGNED, false);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_PEPTIDEMZ, false);
			set = true;
		}
		else if (action->text().toStdString() == "Show by precursor m/z")
		{
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::F_UNASSIGNED, true);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_PEPTIDEMZ, false);
			set = true;
		}
		else if (action->text().toStdString() == "Show by peptide mass")
		{
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::F_UNASSIGNED, true);
      getActive2DWidget()->canvas()->setLayerFlag(LayerData::I_PEPTIDEMZ, true);
			set = true;
		}

		// button is simply pressed
		if (!set)
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
		QAction* action = qobject_cast<QAction *>(sender());
    if (Spectrum2DWidget* win = getActive2DWidget())
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
			//consensus features
			else if (action == dm_elements_2d_)
			{
				 win->canvas()->setLayerFlag(LayerData::C_ELEMENTS,on);
			}
			// identifications
			else if (action == dm_ident_2d_)
			{
				win->canvas()->setLayerFlag(LayerData::I_PEPTIDEMZ, on);
			}
		}
  }

  void TOPPViewBase::updateToolBar()
  {
    SpectrumWidget* w = getActiveSpectrumWidget();

    if (w)
    {
      //set intensity mode
      if(intensity_button_group_->button(w->canvas()->getIntensityMode()))
      {
        intensity_button_group_->button(w->canvas()->getIntensityMode())->setChecked(true);
      } else
      {
        showLogMessage_(LS_ERROR, __PRETTY_FUNCTION__ ,"Button for intensity mode doesn't exist");
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
      if(w2->canvas()->getLayerCount() > 0)
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
          dm_label_2d_->setChecked(w2->canvas()->getCurrentLayer().label!=LayerData::L_NONE);
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
  	//reset
		layer_manager_->clear();
    SpectrumCanvas* cc = getActiveCanvas();
    if (cc == 0) return;

		//determine if this is a 1D view (for text color)
		bool is_1d_view = false;
		if (dynamic_cast<Spectrum1DCanvas*>(cc)) is_1d_view = true;

		layer_manager_->blockSignals(true);
		QListWidgetItem* item = 0;
		QString name;
    for (Size i = 0; i<cc->getLayerCount(); ++i)
    {
    	const LayerData& layer = cc->getLayer(i);
    	//add item
    	item = new QListWidgetItem( layer_manager_ );
			name = layer.name.toQString();
			if (layer.flipped)
			{
				name += " [flipped]";
			}
			item->setText(name);
			if (is_1d_view && cc->getLayerCount()>1)
			{
				QPixmap icon(7,7);
				icon.fill(QColor(layer.param.getValue("peak_color").toQString()));
				item->setIcon(icon);
			}
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

  void TOPPViewBase::updateViewBar()
  {
    SpectrumCanvas* cc = getActiveCanvas();
    int layer_row = layer_manager_->currentRow();

    if (layer_row == -1 || cc == 0)
    {
      // TODO: we need to clean up the SpectraViewWidget & friends
      return;
    }

    if (spectra_view_widget_->isVisible())
    {
      spectra_view_widget_->updateEntries(cc->getCurrentLayer());
    }

    if (spectra_identification_view_widget_->isVisible())
    {
      spectra_identification_view_widget_->attachLayer(&cc->getCurrentLayer());
      spectra_identification_view_widget_->updateEntries();
    }
  }

  void TOPPViewBase::viewChanged(int tab_index)
  {
    // notify that behavior will be deactivated
    view_behavior_->deactivateBehavior();

    // set new behavior
    if (views_tabwidget_->tabText(tab_index) == "Scan view")
    {
      layer_dock_widget_->show();
      filter_dock_widget_->show();
      view_behavior_ = spectraview_behavior_;
    } else if (views_tabwidget_->tabText(tab_index) == "Identification view")
    {
      layer_dock_widget_->show();
      filter_dock_widget_->show();
      if (getActive2DWidget())  // currently 2D window is open
      {
        showSpectrumAs1D(0);
      }
      view_behavior_ = identificationview_behavior_;
    } else if (views_tabwidget_->tabText(tab_index) == "TOPPAS view")
    {
      layer_dock_widget_->hide();
      filter_dock_widget_->hide();
      // no complex behavior for TOPPAS view needed
    } else
    {
      cerr << "Error: tab_index " << tab_index << endl;
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

    // notify that new behavior has been activated
    view_behavior_->activateBehavior();
    updateViewBar();
  }

  /*
  void TOPPViewBase::updateDataBar()
  {
    const set<String> filenames = getFilenamesOfOpenFiles();
    //reset
    data_manager_view_->clear();
    data_manager_view_->blockSignals(true);
    QTreeWidgetItem* item = 0;
    for (set<String>::const_iterator it = filenames.begin(); it != filenames.end(); ++it)
    {
      item = new QTreeWidgetItem(data_manager_view_);
      QString name = it->toQString();
      QFileInfo fi(name);
      item->setText(0, fi.fileName());
    }
    layer_manager_->blockSignals(false);
  }
*/
	void TOPPViewBase::layerSelectionChange(int i)
	{
		if (i!=-1)
		{
      getActiveCanvas()->activateLayer(i);
			updateFilterBar();
      updateViewBar();
		}
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

      if (getActiveCanvas()->getLayer(layer).flipped)
			{
				new_action = context_menu->addAction("Flip upwards (1D)");
			}
			else
			{
				new_action = context_menu->addAction("Flip downwards (1D)");
			}
      if (!getActive1DWidget())
			{
				new_action->setEnabled(false);
			}

			context_menu->addSeparator();
			context_menu->addAction("Preferences");

			QAction* selected = context_menu->exec(layer_manager_->mapToGlobal(pos));
			//delete layer
			if (selected!=0 && selected->text()=="Delete")
			{
        getActiveCanvas()->removeLayer(layer);
			}
			//rename layer
			else if (selected!=0 && selected->text()=="Rename")
			{
				QString name = QInputDialog::getText(this,"Rename layer","Name:");
				if (name!="")
				{
          getActiveCanvas()->setLayerName(layer, name);
				}
			}
			// flip layer up/downwards
			else if (selected != 0 && selected->text() == "Flip downwards (1D)")
			{
        getActive1DWidget()->canvas()->flipLayer(layer);
        getActive1DWidget()->canvas()->setMirrorModeActive(true);
			}
			else if (selected != 0 && selected->text() == "Flip upwards (1D)")
			{
        getActive1DWidget()->canvas()->flipLayer(layer);
        bool b = getActive1DWidget()->canvas()->flippedLayersExist();
        getActive1DWidget()->canvas()->setMirrorModeActive(b);
			}
			else if (selected != 0 && selected->text() == "Preferences")
			{
        getActiveCanvas()->showCurrentLayerPreferences();
			}

			//Update tab bar and window title
      if (getActiveCanvas()->getLayerCount()!=0)
			{
        tab_bar_->setTabText(tab_bar_->currentIndex(), getActiveCanvas()->getLayer(0).name.toQString());
        getActiveSpectrumWidget()->setWindowTitle(getActiveCanvas()->getLayer(0).name.toQString());
			}
			else
			{
				tab_bar_->setTabText(tab_bar_->currentIndex(),"empty");
        getActiveSpectrumWidget()->setWindowTitle("empty");
			}

			//Update filter bar, spectrum bar and layer bar
			updateLayerBar();
      updateViewBar();
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
    if (selected != 0 && selected->text()== "Clear")
		{
			log_->clear();
		}
		delete (context_menu);
	}

	void TOPPViewBase::filterContextMenu(const QPoint & pos)
	{
		//do nothing if no window is open
    if (getActiveCanvas()==0) return;

    //do nothing if no layer is loaded into the canvas
    if (getActiveCanvas()->getLayerCount() == 0) return;

		QMenu* context_menu = new QMenu(filters_);

		//warn if the current layer is not visible
    String layer_name = String("Layer: ") + getActiveCanvas()->getCurrentLayer().name;
    if (!getActiveCanvas()->getCurrentLayer().visible)
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
        DataFilters filters = getActiveCanvas()->getCurrentLayer().filters;
				filters.remove(filters_->row(item));
        getActiveCanvas()->setFilters(filters);
				updateFilterBar();
			}
			else if (selected->text()=="Edit")
			{
				filterEdit(item);
			}
			else if (selected->text()=="Add filter")
			{
        DataFilters filters = getActiveCanvas()->getCurrentLayer().filters;
				DataFilters::DataFilter filter;
				DataFilterDialog dlg(filter, this);
				if (dlg.exec())
				{
					filters.add(filter);
          getActiveCanvas()->setFilters(filters);
					updateFilterBar();
				}
			}
		}
		delete (context_menu);
	}

	void TOPPViewBase::filterEdit(QListWidgetItem* item)
	{
    DataFilters filters = getActiveCanvas()->getCurrentLayer().filters;
		DataFilters::DataFilter filter = filters[filters_->row(item)];
		DataFilterDialog dlg(filter, this);
		if (dlg.exec())
		{
			filters.replace(filters_->row(item),filter);
      getActiveCanvas()->setFilters(filters);
			updateFilterBar();
		}
	}

	void TOPPViewBase::layerEdit(QListWidgetItem* /*item*/)
	{
    getActiveCanvas()->showCurrentLayerPreferences();
	}

  void TOPPViewBase::updateFilterBar()
  {
  	//update filters
  	filters_->clear();

    SpectrumCanvas* canvas = getActiveCanvas();
		if (canvas==0) return;
		if (canvas->getLayerCount()==0) return;

    const DataFilters& filters = getActiveCanvas()->getCurrentLayer().filters;
		for (Size i=0; i<filters.size(); ++i)
		{
			QListWidgetItem* item = new QListWidgetItem(filters_);
			item->setText(filters[i].toString().toQString());
		}

  	//update check box
    filters_check_box_->setChecked(getActiveCanvas()->getCurrentLayer().filters.isActive());
  }

	void TOPPViewBase::layerFilterVisibilityChange(bool on)
	{
    if (getActiveCanvas())
		{
      getActiveCanvas()->changeLayerFilterState(getActiveCanvas()->activeLayerIndex(),on);
		}
	}

	void TOPPViewBase::layerVisibilityChange(QListWidgetItem* item)
	{
		int layer;
		bool visible;
		layer = layer_manager_->row(item);
    visible = getActiveCanvas()->getLayer(layer).visible;

		if (item->checkState()==Qt::Unchecked && visible)
		{
      getActiveCanvas()->changeVisibility(layer, false);
		}
		else if (item->checkState()==Qt::Checked && !visible)
		{
      getActiveCanvas()->changeVisibility(layer, true);
		}
	}

  void TOPPViewBase::updateTabBar(QWidget* w)
  {
  	if (w)
  	{
      EnhancedTabBarWidgetInterface* tbw = dynamic_cast<EnhancedTabBarWidgetInterface*>(w);
      Int window_id = tbw->getWindowId();
  		tab_bar_->setCurrentId(window_id);
  	}
  }

  void TOPPViewBase::tileVertical()
  {
    // primitive horizontal tiling
    QWidgetList windows = ws_->windowList();
    if ( !windows.count() ) return;

    if (getActive1DWidget()) getActive1DWidget()->showNormal();
    if (getActive2DWidget()) getActive2DWidget()->showNormal();

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

    if (getActive1DWidget()) getActive1DWidget()->showNormal();
    if (getActive2DWidget()) getActive2DWidget()->showNormal();

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

  void TOPPViewBase::linkZoom()
  {
    zoom_together_ = !zoom_together_;
    if(!zoom_together_)
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
    updateToolBar();
    updateViewBar();
    updateCurrentPath();
  }

  void TOPPViewBase::layerZoomChanged()
  {
    QWidgetList windows = ws_->windowList();
    if ( !windows.count() ) return;
    if ( !zoom_together_ ) return;

    SpectrumWidget* w = getActiveSpectrumWidget();

    // figure out which dimension the active widget has: 2D (MSExperiment) or 1D (Iontrace)
    // and get the corresponding RT values.
    Spectrum1DWidget* sw1 = qobject_cast<Spectrum1DWidget*>(w);
    Spectrum2DWidget* sw2 = qobject_cast<Spectrum2DWidget*>(w);
    Spectrum3DWidget* sw3 = qobject_cast<Spectrum3DWidget*>(w);
    int widget_dimension = -1;
    if (sw1 != 0)
    {
      widget_dimension = 1;
    }
    else if (sw2 != 0)
    {
      widget_dimension = 2;
    }
    else if (sw3 != 0)
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
    if(getActiveCanvas()->getCurrentLayer().type == LayerData::DT_CHROMATOGRAM || 
        (getActiveCanvas()->getCurrentLayer().getPeakData()->size() > 0 && 
         getActiveCanvas()->getCurrentLayer().getPeakData()->metaValueExists("is_chromatogram") && 
         getActiveCanvas()->getCurrentLayer().getPeakData()->getMetaValue("is_chromatogram").toBool() 
         ) ) {
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
      for ( int i = 0; i < int(windows.count()); ++i )
      {
        QWidget *window = windows.at(i);
        DRange<2> visible_area;
        SpectrumWidget* specwidg = qobject_cast<SpectrumWidget*>(window);

        // Skip if its not a SpectrumWidget, if it is not a chromatogram or if the dimensions don't match.
        if(!specwidg) continue;
        if(!(specwidg->canvas()->getCurrentLayer().type == LayerData::DT_CHROMATOGRAM) && 
           !(specwidg->canvas()->getCurrentLayer().getPeakData()->size() > 0 && 
             specwidg->canvas()->getCurrentLayer().getPeakData()->metaValueExists("is_chromatogram") && 
             specwidg->canvas()->getCurrentLayer().getPeakData()->getMetaValue("is_chromatogram").toBool() 
             ) ) 
        {
          continue;
        }
        if( !(widget_dimension == 1 && qobject_cast<Spectrum1DWidget*>(specwidg)) && 
            !(widget_dimension == 2 && qobject_cast<Spectrum2DWidget*>(specwidg))  )
        {
          continue;
        }

        visible_area = specwidg->canvas()->getVisibleArea();

        // if we found a min/max RT, change all windows of 1 dimension
        if(minRT != -1 && maxRT != -1 && qobject_cast<Spectrum1DWidget*>(window))
        {
          visible_area.setMinX(minRT);
          visible_area.setMaxX(maxRT);
        }
        specwidg->canvas()->setVisibleArea(visible_area);
      }
    }
    else {
      DRange<2> new_visible_area = w->canvas()->getVisibleArea();
      // go through all windows, adjust the visible area where necessary
      for(int i = 0; i < int(windows.count()); ++i)
      {
        QWidget *window = windows.at(i);
        SpectrumWidget* specwidg = qobject_cast<SpectrumWidget*>(window);

        // Skip if its not a SpectrumWidget, if it is a chromatogram or if the dimensions don't match.
        if(!specwidg) continue;
        if((specwidg->canvas()->getCurrentLayer().type == LayerData::DT_CHROMATOGRAM) || 
           (specwidg->canvas()->getCurrentLayer().getPeakData()->size() > 0 && 
             specwidg->canvas()->getCurrentLayer().getPeakData()->metaValueExists("is_chromatogram") && 
             specwidg->canvas()->getCurrentLayer().getPeakData()->getMetaValue("is_chromatogram").toBool() 
             ) )
        {
          continue;
        }
        if( !(widget_dimension == 1 && qobject_cast<Spectrum1DWidget*>(specwidg)) && 
            !(widget_dimension == 2 && qobject_cast<Spectrum2DWidget*>(specwidg))  )
        {
          continue;
        }
        specwidg->canvas()->setVisibleArea(new_visible_area);
      }
      return;
    }

  }

  void TOPPViewBase::layerDeactivated()
  {

  }

  void TOPPViewBase::showSpectrumWidgetInWindow(SpectrumWidget* sw, const String& caption)
  {
  	ws_->addWindow(sw);
    connect(sw->canvas(),SIGNAL(preferencesChange()),this,SLOT(updateLayerBar()));
    connect(sw->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(layerActivated()));
    connect(sw->canvas(),SIGNAL(layerModficationChange(Size,bool)),this,SLOT(updateLayerBar()));
    connect(sw->canvas(),SIGNAL(layerZoomChanged(QWidget*)),this,SLOT(layerZoomChanged()));
    connect(sw,SIGNAL(sendStatusMessage(std::string,OpenMS::UInt)),this,SLOT(showStatusMessage(std::string,OpenMS::UInt)));
    connect(sw,SIGNAL(sendCursorStatus(double,double)),this,SLOT(showCursorStatus(double,double)));
    connect(sw,SIGNAL(dropReceived(const QMimeData*,QWidget*,int)),this,SLOT(copyLayer(const QMimeData*,QWidget*,int)));

    // 1D spectrum specific signals
    Spectrum1DWidget* sw1 = qobject_cast<Spectrum1DWidget*>(sw);
    if (sw1 != 0)
    {
      connect(sw1, SIGNAL(showCurrentPeaksAs2D()), this, SLOT(showCurrentPeaksAs2D()));
      connect(sw1, SIGNAL(showCurrentPeaksAs3D()), this, SLOT(showCurrentPeaksAs3D()));
    }

    // 2D spectrum specific signals
  	Spectrum2DWidget* sw2 = qobject_cast<Spectrum2DWidget*>(sw);
  	if (sw2 != 0)
  	{
  		connect(sw2->getHorizontalProjection(),SIGNAL(sendCursorStatus(double,double)),this,SLOT(showCursorStatus(double,double)));
  		connect(sw2->getVerticalProjection(),SIGNAL(sendCursorStatus(double,double)),this,SLOT(showCursorStatusInvert(double,double)));
      connect(sw2, SIGNAL(showSpectrumAs1D(int)), this, SLOT(showSpectrumAs1D(int)));
      connect(sw2, SIGNAL(showSpectrumAs1D(std::vector<int, std::allocator<int> >)), this, SLOT(showSpectrumAs1D(std::vector<int, std::allocator<int> >)));
      connect(sw2, SIGNAL(showCurrentPeaksAs3D()), this, SLOT(showCurrentPeaksAs3D()));
  	}

    // 3D spectrum specific signals
    Spectrum3DWidget* sw3 = qobject_cast<Spectrum3DWidget*>(sw);
    if (sw3 != 0)
    {
      connect(sw3, SIGNAL(showCurrentPeaksAs2D()), this, SLOT(showCurrentPeaksAs2D()));
    }

	  sw->setWindowTitle(caption.toQString());

		//add tab with id
  	static int window_counter = 4711;

    sw->setWindowId(window_counter++);

    tab_bar_->addTab(caption.toQString(), sw->getWindowId());

    //connect slots and sigals for removing the widget from the bar, when it is closed
    //- through the menu entry
    //- through the tab bar
    //- thourgh the MDI close button
    connect(sw,SIGNAL(aboutToBeDestroyed(int)),tab_bar_,SLOT(removeId(int)));

    tab_bar_->setCurrentId(sw->getWindowId());

		//show first window maximized (only visible windows are in the list)
		if (ws_->windowList().count()==0)
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

  void TOPPViewBase::activate1DSpectrum(int index)
  {
    Spectrum1DWidget* w = getActive1DWidget();
    if (w)
    {
      view_behavior_->activate1DSpectrum(index);
    }
  }
  void TOPPViewBase::activate1DSpectrum(std::vector<int, std::allocator<int> > indices)
  {
    Spectrum1DWidget* w = getActive1DWidget();
    if (w)
    {
      view_behavior_->activate1DSpectrum(indices);
    }
  }

  void TOPPViewBase::deactivate1DSpectrum(int index)
  {
    Spectrum1DWidget* w = getActive1DWidget();
    if (w)
    {
      view_behavior_->deactivate1DSpectrum(index);
    }
  }

  EnhancedWorkspace* TOPPViewBase::getWorkspace() const
  {
    return ws_;
  }

  SpectrumWidget* TOPPViewBase::getActiveSpectrumWidget() const
  {
    if (!ws_->activeWindow())
    {
      return 0;
    }
    return qobject_cast<SpectrumWidget*>(ws_->activeWindow());
  }

  SpectrumCanvas*  TOPPViewBase::getActiveCanvas() const
  {
    SpectrumWidget* sw = qobject_cast<SpectrumWidget*>(ws_->activeWindow());
    if (sw == 0)
    {
    	return 0;
    }
    return sw->canvas();
  }


  TOPPASWidget* TOPPViewBase::getActiveTOPPASWidget() const
  {
    if (!ws_->activeWindow())
    {
      return 0;
    }
    return qobject_cast<TOPPASWidget*>(ws_->activeWindow());
  }

  Spectrum1DWidget* TOPPViewBase::getActive1DWidget() const
  {
    Spectrum1DWidget* w = qobject_cast<Spectrum1DWidget*>(getActiveSpectrumWidget());
    if (!w)
    {
      return 0;
    }
		return w;
  }

  Spectrum2DWidget* TOPPViewBase::getActive2DWidget() const
  {
    Spectrum2DWidget* w = qobject_cast<Spectrum2DWidget*>(getActiveSpectrumWidget());
    if (!w)
    {
      return 0;
    }
		return w;
  }

  Spectrum3DWidget* TOPPViewBase::getActive3DWidget() const
  {
    Spectrum3DWidget* w = qobject_cast<Spectrum3DWidget*>(getActiveSpectrumWidget());
    if (!w)
    {
      return 0;
    }
		return w;
  }

  void TOPPViewBase::loadPreferences(String filename)
  {
    // compose default ini file path
    String default_ini_file = String(QDir::homePath()) + "/.TOPPView.ini";

    if (filename == "")
    {
      filename = default_ini_file;
    }

    // load preferences, if file exists
    if (File::exists(filename))
    {
    	bool error = false;
    	Param tmp;
      try
      { // the file might be corrupt
    	  tmp.load(filename);
      }
      catch (...)
      {
        error = true;
      }

      //apply preferences if they are of the current TOPPView version
      if(!error && tmp.exists("preferences:version") &&
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
    catch(Exception::UnableToCreateFile& /*e*/)
    {
      cerr << "Unable to create INI File: '" << string(param_.getValue("PreferencesFile")) << "'" << endl;
    }
  }

  void TOPPViewBase::openRecentFile()
  {
		QAction* action = qobject_cast<QAction *>(sender());
    if (action)
		{
      QString filename = action->text();
      if (filename.endsWith(".toppas", Qt::CaseInsensitive))
      {
        addTOPPASFile(filename, true);
      }
      else
      {
        addDataFile(filename, true, true);
      }
    }
	}

  QStringList TOPPViewBase::getFileList_(const String& path_overwrite)
  {
    String filter_all = "readable files (*.mzML *.mzXML *.mzData *.featureXML *.consensusXML *.idXML *.dta *.dta2d fid *.bz2 *.gz *.toppas);;";
    String filter_single = "mzML files (*.mzML);;mzXML files (*.mzXML);;mzData files (*.mzData);;feature map (*.featureXML);;consensus feature map (*.consensusXML);;peptide identifications (*.idXML);;XML files (*.xml);;XMass Analysis (fid);;dta files (*.dta);;dta2d files (*.dta2d);;bzipped files (*.bz2);;gzipped files (*.gz);;TOPPAS files (*.toppas);;all files (*)";

		QString open_path = current_path_.toQString();
		if (path_overwrite!="")
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

		return file_names;
  }

  void TOPPViewBase::openFileDialog()
  {
	 	QStringList files = getFileList_();
    for(QStringList::iterator it = files.begin(); it!=files.end(); ++it)
		{
      QString filename = *it;
      if (filename.endsWith(".toppas", Qt::CaseInsensitive))
      {
        addTOPPASFile(filename, true);
      }
      else
      {
        addDataFile(filename, true, true);
      }
		}
  }

  void TOPPViewBase::addTOPPASFile(const String& file_name, bool in_new_window)
  {
    TOPPASScene* scene = 0;
    if (in_new_window) // open in new window
    {
      // create TOPPASWidget, load data and open in new window
      TOPPASWidget* tw = new TOPPASWidget(Param(), ws_, toppas_tmp_path_);
      showTOPPipelineInWindow_(tw, File::basename(file_name));
      scene = tw->getScene();
      scene->load(file_name);
      addRecentFile_(file_name);
    }
    else  // merge into existing scene
    {
      TOPPASWidget* tw = getActiveTOPPASWidget();
      if (!tw)
      {
        return;
      }
      // create TOPPASWidget, load data into temporary scene and include into existing one
      TOPPASScene* tmp_scene = new TOPPASScene(0, File::getTempDirectory().toQString()+QDir::separator(), false);
      tmp_scene->load(file_name);
      scene = tw->getScene();
      scene->include(tmp_scene);
      delete tmp_scene;
    }

    if (in_new_window)
    {
      // connect scene signals only if we created a new window (otherwise they already exist)
      connect(scene, SIGNAL(saveMe()), this, SLOT(savePipeline()));
      connect(scene, SIGNAL(selectionCopied(TOPPASScene*)), this, SLOT(saveToClipboard(TOPPASScene*)));
      connect(scene, SIGNAL(requestClipboardContent()), this, SLOT(sendClipboardContent()));
      connect(scene, SIGNAL(mainWindowNeedsUpdate()), this, SLOT(updateMenu()));
      connect(scene, SIGNAL(openInTOPPView(QStringList)), this, SLOT(openFilesInTOPPView(QStringList)));
    }

    //connect vertex signals/slots for log messages
    for (TOPPASScene::VertexIterator it = scene->verticesBegin(); it != scene->verticesEnd(); ++it)
    {
      TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(*it);
      if (tv)
      {
        connect(tv, SIGNAL(toolStarted()), this, SLOT(toolStarted()));
        connect(tv, SIGNAL(toolFinished()), this, SLOT(toolFinished()));
        connect(tv, SIGNAL(toolCrashed()), this, SLOT(toolCrashed()));
        connect(tv, SIGNAL(toolFailed()), this, SLOT(toolFailed()));
        connect(tv, SIGNAL(toppOutputReady(const QString&)), this, SLOT(updateTOPPOutputLog(const QString&)));
        continue;
      }
      TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(*it);
      if (oflv)
      {
        connect (oflv, SIGNAL(outputFileWritten(const String&)), this, SLOT(outputVertexFinished(const String&)));
        continue;
      }
    }
  }

  void TOPPViewBase::showTOPPipelineInWindow_(TOPPASWidget* tw, const String& caption)
  {
    ws_->addWindow(tw);

    connect(tw,SIGNAL(sendStatusMessage(std::string,OpenMS::UInt)),this,SLOT(showStatusMessage(std::string,OpenMS::UInt)));
    connect(tw,SIGNAL(sendCursorStatus(double,double)),this,SLOT(showCursorStatus(double,double)));
    connect(tw,SIGNAL(toolDroppedOnWidget(double,double)),this,SLOT(insertNewVertex_(double,double)));
    connect(tw,SIGNAL(pipelineDroppedOnWidget(const String&, bool)),this,SLOT(addTOPPASFile(const String&, bool)));

    tw->setWindowTitle(caption.toQString());

    //add tab with id
    static int window_counter = 1337;
    tw->setWindowId(window_counter);
    window_counter++;

    //connect slots and signals for removing the widget from the bar, when it is closed
    //- through the menu entry
    //- through the tab bar
    //- through the MDI close button
    connect(tw, SIGNAL(aboutToBeDestroyed(int)),tab_bar_,SLOT(removeId(int)));

    //show first window maximized (only visible windows are in the list)
    if (ws_->windowList().count()==0)
    {
      tw->showMaximized();
    }
    else
    {
      tw->show();
    }
    TOPPASScene* ts = tw->getScene();

    connect (ts, SIGNAL(entirePipelineFinished()), this, SLOT(showPipelineFinishedLogMessage()));
    connect (ts, SIGNAL(entirePipelineFinished()), this, SLOT(updateMenu()));
    connect (ts, SIGNAL(pipelineExecutionFailed()), this, SLOT(updateMenu()));

    ts->setSceneRect((tw->mapToScene(tw->rect())).boundingRect());

    tab_bar_->addTab(caption.toQString(), tw->getWindowId());
    tab_bar_->setCurrentId(tw->getWindowId());
    enhancedWorkspaceWindowChanged(tw->getWindowId());
  }

  void TOPPViewBase::insertNewVertex_(double x, double y, QTreeWidgetItem* item)
  {
    // get toppas tree view from tab widget 2
    TOPPASTreeView* toppas_tree_view = qobject_cast<TOPPASTreeView*>(views_tabwidget_->widget(2));

    if (!getActiveTOPPASWidget() || !getActiveTOPPASWidget()->getScene() || !toppas_tree_view)
    {
      return;
    }

    TOPPASScene* scene = getActiveTOPPASWidget()->getScene();
    QTreeWidgetItem* current_tool = item ? item : toppas_tree_view->currentItem();
    String tool_name = String(current_tool->text(0));
    TOPPASVertex* tv = 0;

    if (tool_name == "<Input files>")
    {
      tv = new TOPPASInputFileListVertex();
    }
    else if (tool_name == "<Output files>")
    {
      tv = new TOPPASOutputFileListVertex();
      TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(tv);
      connect (oflv, SIGNAL(outputFileWritten(const String&)), this, SLOT(outputVertexFinished(const String&)));
      scene->connectOutputVertexSignals(oflv);
    }
    else if (tool_name == "<Merger>")
    {
      tv = new TOPPASMergerVertex(true);
    }
    else if (tool_name == "<Collector>")
    {
      tv = new TOPPASMergerVertex(false);
    }
    else // node is a TOPP tool
    {
      if (current_tool->childCount() > 0)
      {
        // category or tool name with types is selected (instead of a concrete type)
        return;
      }
      String tool_type;
      if (current_tool->parent() != 0 && current_tool->parent()->parent() != 0)
      {
        // selected item is a type
        tool_type = String(current_tool->text(0));
        tool_name = String(current_tool->parent()->text(0));
      }
      else
      {
        // normal tool which does not have type selected
        tool_name = String(current_tool->text(0));
        tool_type = "";
      }

      tv = new TOPPASToolVertex(tool_name, tool_type);
      TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(tv);
      connect (ttv, SIGNAL(toolStarted()), this, SLOT(toolStarted()));
      connect (ttv, SIGNAL(toolFinished()), this, SLOT(toolFinished()));
      connect (ttv, SIGNAL(toolCrashed()), this, SLOT(toolCrashed()));
      connect (ttv, SIGNAL(toolFailed()), this, SLOT(toolFailed()));
      connect (ttv, SIGNAL(toppOutputReady(const QString&)), this, SLOT(updateTOPPOutputLog(const QString&)));

      scene->connectToolVertexSignals(ttv);
    }

    scene->connectVertexSignals(tv);
    scene->addVertex(tv);
    tv->setPos(x,y);
    tv->setZValue(toppas_z_value_);
    toppas_z_value_ += 0.000001;
    scene->topoSort();
    scene->setChanged(true);
  }

  void TOPPViewBase::insertNewVertexInCenter_(QTreeWidgetItem* item)
  {
    TOPPASTreeView* toppas_tree_view = qobject_cast<TOPPASTreeView*>(views_tabwidget_->widget(2));

    if (!getActiveTOPPASWidget() || !getActiveTOPPASWidget()->getScene() || !toppas_tree_view || !toppas_tree_view->currentItem())
    {
      return;
    }

    QPointF insert_pos = getActiveTOPPASWidget()->mapToScene(QPoint((getActiveTOPPASWidget()->width()/2.0)+(qreal)(5*toppas_node_offset_),(getActiveTOPPASWidget()->height()/2.0)+(qreal)(5*toppas_node_offset_)));
    insertNewVertex_(insert_pos.x(), insert_pos.y(), item);
    toppas_node_offset_ = (toppas_node_offset_+1) % 10;
  }

  void TOPPViewBase::openExampleDialog()
  {
	 	QStringList files = getFileList_(File::getOpenMSDataPath() + "/examples/");

    for(QStringList::iterator it = files.begin(); it != files.end(); ++it)
		{
      QString filename = *it;
      if (filename.endsWith(".toppas", Qt::CaseInsensitive))
      {
        addTOPPASFile(filename, true);
      }
      else
      {
        addDataFile(filename, true, true);
      }
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
			catch (DBConnection::InvalidQuery& er)
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
          addDataDB(*it, true);
				}
			}
		}
  }

	void TOPPViewBase::rerunTOPPTool()
	{
		//warn if hidden layer => wrong layer selected...
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();
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
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();
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
		ToolsDialog tools_dialog(this,topp_.file_name+"_ini",current_path_,getCurrentLayer()->type);

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
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();

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
				f.store(topp_.file_name+"_in",exp);
			}
			else
			{
        f.store(topp_.file_name+"_in",*layer.getPeakData());
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
        f.store(topp_.file_name+"_in",exp);
      }
      else
      {
        f.store(topp_.file_name+"_in",*layer.getPeakData());
      }
    }
    else if (layer.type==LayerData::DT_FEATURE)
		{
			if (topp_.visible)
			{
				FeatureMapType map;
        getActiveCanvas()->getVisibleFeatureData(map);
				FeatureXMLFile().store(topp_.file_name+"_in",map);
			}
			else
			{
        FeatureXMLFile().store(topp_.file_name+"_in",*layer.getFeatureMap());
			}
		}
		else
		{
			if (topp_.visible)
			{
				ConsensusMapType map;
        getActiveCanvas()->getVisibleConsensusData(map);
				ConsensusXMLFile().store(topp_.file_name+"_in",map);
			}
			else
			{
        ConsensusXMLFile().store(topp_.file_name+"_in",*layer.getConsensusMap());
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
				addDataFile(topp_.file_name+"_out",true,false, topp_.layer_name + " (" + topp_.tool + ")", topp_.window_id, topp_.spectrum_id);
			}
		}

		//clean up
		delete topp_.process;
		topp_.process = 0;
		updateMenu();

		//clean up temporary files
		if (param_.getValue("preferences:topp_cleanup")=="true")
		{
			File::remove(topp_.file_name+"_ini");
			File::remove(topp_.file_name+"_in");
			File::remove(topp_.file_name+"_out");
  	}
  }

	const LayerData* TOPPViewBase::getCurrentLayer() const
	{
    SpectrumCanvas* canvas = getActiveCanvas();
		if (canvas==0)
		{
			return 0;
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
    		setMinimumSize(700,700);
    	}
    	else
    	{
    		setMinimumSize(400,400);
    	}
    	w->toggleProjections();
    }
  }

	void TOPPViewBase::annotateWithID()
	{
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();
		//warn if hidden layer => wrong layer selected...
		if (!layer.visible)
		{
  		showLogMessage_(LS_NOTICE,"The current layer is not visible","Have you selected the right layer for this action?");
		}

		//load id data
    QString name = QFileDialog::getOpenFileName(this,
                                                "Select protein identification data",
                                                current_path_.toQString(),
                                                "idXML files (*.idXML);; all files (*.*)");

		if(name!="")
		{
			vector<PeptideIdentification> identifications;
			vector<ProteinIdentification> protein_identifications;
  
      try
      {
  			String document_id;
  			IdXMLFile().load(name, protein_identifications, identifications, document_id);
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::warning(this, "Error", QString("Loading of idXML file failed! (") + e.what() + ")");
			  return;
      }

      IDMapper mapper;
			if (layer.type==LayerData::DT_PEAK)
			{
        Param p = mapper.getDefaults();
        p.setValue("rt_tolerance", 0.1, "RT tolerance (in seconds) for the matching");
        p.setValue("mz_tolerance", 1.0, "m/z tolerance (in ppm or Da) for the matching");
        p.setValue("mz_measure", "Da", "unit of 'mz_tolerance' (ppm or Da)");
        mapper.setParameters(p);
        mapper.annotate(*layer.getPeakData(), identifications, protein_identifications);
        views_tabwidget_->setTabEnabled(1, true);         // enable identification view
			}
			else if (layer.type==LayerData::DT_FEATURE)
			{
        mapper.annotate(*layer.getFeatureMap(), identifications, protein_identifications);
			}
			else
			{
        mapper.annotate(*layer.getConsensusMap(), identifications, protein_identifications);
			}
		}
    updateViewBar();
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
      AASequence aa_sequence;
      try
      {
        aa_sequence.setStringSequence(seq_string);
      }
      catch (Exception::BaseException& e)
      {
        QMessageBox::warning(this, "Error", QString("Spectrum generation failed! (") + e.what() + ")");
			  return;
      }
			Int charge = spec_gen_dialog.spin_box->value();

			if (aa_sequence.isValid())
			{
				RichPeakSpectrum rich_spec;
				TheoreticalSpectrumGenerator generator;
				Param p;

        p.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");

				bool losses = (spec_gen_dialog.list_widget->item(7)->checkState() == Qt::Checked); // "Neutral losses"
				String losses_str = losses ? "true" : "false";
				p.setValue("add_losses", losses_str, "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");

				bool isotopes = (spec_gen_dialog.list_widget->item(8)->checkState() == Qt::Checked); // "Isotope clusters"
				String iso_str = isotopes ? "true" : "false";
				p.setValue("add_isotopes", iso_str, "If set to 1 isotope peaks of the product ion peaks are added");

        bool abundant_immonium_ions = (spec_gen_dialog.list_widget->item(9)->checkState() == Qt::Checked); // "abundant immonium-ions"
        String abundant_immonium_ions_str = abundant_immonium_ions ? "true" : "false";
        p.setValue("add_abundant_immonium_ions", abundant_immonium_ions_str, "Add most abundant immonium ions");

        Size max_iso_count = (Size)spec_gen_dialog.max_iso_spinbox->value();
        p.setValue("max_isotope", max_iso_count, "Number of isotopic peaks");
				p.setValue("a_intensity", spec_gen_dialog.a_intensity->value(), "Intensity of the a-ions");
				p.setValue("b_intensity", spec_gen_dialog.b_intensity->value(), "Intensity of the b-ions");
				p.setValue("c_intensity", spec_gen_dialog.c_intensity->value(), "Intensity of the c-ions");
				p.setValue("x_intensity", spec_gen_dialog.x_intensity->value(), "Intensity of the x-ions");
				p.setValue("y_intensity", spec_gen_dialog.y_intensity->value(), "Intensity of the y-ions");
				p.setValue("z_intensity", spec_gen_dialog.z_intensity->value(), "Intensity of the z-ions");
				DoubleReal rel_loss_int = (DoubleReal)(spec_gen_dialog.rel_loss_intensity->value()) / 100.0;
				p.setValue("relative_loss_intensity", rel_loss_int, "Intensity of loss ions, in relation to the intact ion intensity");
				generator.setParameters(p);

        try
        {
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
          if (spec_gen_dialog.list_widget->item(9)->checkState() == Qt::Checked) // "abundant Immonium-ions"
          {
            generator.addAbundantImmoniumIons(rich_spec);
          }
        }
        catch (Exception::BaseException& e)
        {
          QMessageBox::warning(this, "Error", QString("Spectrum generation failed! (") + e.what() + "). Please report this to the developers (specify what input you used)!");
				  return;
        }

        // convert rich spectrum to simple spectrum
				PeakSpectrum new_spec;
				for (RichPeakSpectrum::Iterator it = rich_spec.begin(); it != rich_spec.end(); ++it)
				{
					new_spec.push_back(static_cast<Peak1D>(*it));
				}

				PeakMap new_exp;
				new_exp.push_back(new_spec);
        ExperimentSharedPtrType new_exp_sptr(new PeakMap(new_exp));
        FeatureMapSharedPtrType f_dummy(new FeatureMapType());
        ConsensusMapSharedPtrType c_dummy(new ConsensusMapType());
				vector<PeptideIdentification> p_dummy;
        addData(f_dummy, c_dummy, p_dummy, new_exp_sptr, LayerData::DT_CHROMATOGRAM, false, true, true, "", seq_string + QString(" (theoretical)"));

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
			DoubleReal tolerance = spec_align_dialog.tolerance_spinbox->value();
      param.setValue("tolerance", tolerance, "Defines the absolut (in Da) or relative (in ppm) mass tolerance");
			String unit_is_ppm = spec_align_dialog.ppm->isChecked() ? "true" : "false";
      param.setValue("is_relative_tolerance", unit_is_ppm, "If true, the mass tolerance is interpreted as ppm value otherwise in Dalton");

			active_1d_window->performAlignment((UInt)layer_index_1, (UInt)layer_index_2, param);

			DoubleReal al_score = cc->getAlignmentScore();
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
        spectraview_behavior_->showSpectrumAs1D(index);
      }

      if (spectra_identification_view_widget_->isVisible())
      {
        identificationview_behavior_->showSpectrumAs1D(index);
      }
    } else if (widget_2d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_behavior_->showSpectrumAs1D(index);
      }

      if (spectra_identification_view_widget_->isVisible())
      {
        identificationview_behavior_->showSpectrumAs1D(index);
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
        spectraview_behavior_->showSpectrumAs1D(indices);
      }


    } else if (widget_2d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_behavior_->showSpectrumAs1D(indices);
      }


    }

  }

  void TOPPViewBase::showCurrentPeaksAs2D()
  {
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();
    ExperimentSharedPtrType exp_sptr = layer.getPeakData();

    //open new 2D widget
    Spectrum2DWidget* w = new Spectrum2DWidget(getSpectrumParameters(2), ws_);

    //add data
    if (!w->canvas()->addLayer(exp_sptr, layer.filename))
    {
      return;
    }

    String caption = layer.name;
    w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);
    showSpectrumWidgetInWindow(w, caption);
    updateLayerBar();
    updateViewBar();
    updateFilterBar();
    updateMenu();
  }

  void TOPPViewBase::showCurrentPeaksAs3D()
  {
    const LayerData& layer = getActiveCanvas()->getCurrentLayer();

    if (layer.type == LayerData::DT_PEAK)
    {
      //open new 3D widget
      Spectrum3DWidget* w = new Spectrum3DWidget(getSpectrumParameters(3), ws_);

      ExperimentSharedPtrType exp_sptr = getActiveCanvas()->getCurrentLayer().getPeakData();

      if (!w->canvas()->addLayer(exp_sptr, layer.filename))
      {
        return;
      }

     if (getActive1DWidget()) // switch from 1D to 3D
     {
       //TODO:
       //- doesnt make sense for fragment scan
       //- build new Area with mz range equal to 1D visible range
       //- rt range either overall MS1 data range or some convenient window

     } else if (getActive2DWidget())  // switch from 2D to 3D
     {
        w->canvas()->setVisibleArea(getActiveCanvas()->getVisibleArea());
     }

      // set layer name
      String caption = layer.name + " (3D)";
      w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);
      showSpectrumWidgetInWindow(w, caption);

      // set intensity mode (after spectrum has been added!)
      setIntensityMode(SpectrumCanvas::IM_SNAP);

      updateLayerBar();
      updateViewBar();
      updateFilterBar();
      updateMenu();
    }
    else
    {
      showLogMessage_(LS_NOTICE,"Wrong layer type","You cannot open feature data in 3D mode.");
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
		label->setPixmap(QPixmap(":/TOPP_about.png"));
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

		//close button
		QPushButton* button = new QPushButton("Close",dlg);
		grid->addWidget(button,1,1,Qt::AlignBottom | Qt::AlignRight);
		connect(button,SIGNAL(clicked()),dlg,SLOT(close()));

		//execute
		dlg->exec();
	}

	void TOPPViewBase::updateProcessLog()
	{
		//show log if there is output
		qobject_cast<QWidget *>(log_->parent())->show();

		//update log_
    log_->moveCursor(QTextCursor::End,QTextCursor::MoveAnchor); // move cursor to end, since text is inserted at cursor
    log_->insertPlainText(topp_.process->readAllStandardOutput());
		
	}

  Param TOPPViewBase::getSpectrumParameters(UInt dim)
	{
		Param out = param_.copy(String("preferences:") + dim + "d:",true);
		out.setValue("default_path",param_.getValue("preferences:default_path").toString());
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
    if (getActiveCanvas()!=0)
  	{
  		canvas_exists = true;
  	}
  	//is there a layer?
  	bool layer_exists = false;
    if (canvas_exists && getActiveCanvas()->getLayerCount()!=0)
  	{
  		layer_exists = true;
  	}
		//is there a TOPP tool running
		bool topp_running = false;
		if (topp_.process!=0)
		{
			topp_running = true;
		}

    TOPPASWidget* tw = getActiveTOPPASWidget();
    TOPPASScene* ts = 0;
    if (tw)
    {
      ts = tw->getScene();
    } else // active widget is no TOPPAS widget
    {

    }

    bool mirror_mode = getActive1DWidget() && getActive1DWidget()->canvas()->mirrorModeActive();
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
      else if (text=="&Go to" || text=="&Edit meta data" || text=="&Statistics" || text=="&Annotate with identification"  || text=="Save all data"  || text=="Save visible data"  || text=="Preferences")
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
      else if (text=="&Run (F5)")  // pipeline menu
      {
        bool show = false;
        if (ts && !(ts->isPipelineRunning()))
        {
          show = true;
        }
        actions[i]->setEnabled(show);
      }
      else if (text=="&Abort") // pipeline menu
      {
        bool show = false;
        if (ts && ts->isPipelineRunning())
        {
          show = true;
        }
        actions[i]->setEnabled(show);
      }
      else if (text=="&Include TOPPAS pipeline") // pipeline menu
      {
        bool show = ts;
        actions[i]->setEnabled(show);
      }
      else if (text=="&Load TOPPAS resource file") // pipeline menu
      {
        bool show = ts;
        actions[i]->setEnabled(show);
      }
      else if (text=="Save TOPPAS &resource file") // pipeline menu
      {
        bool show = ts;
        actions[i]->setEnabled(show);
      }
      else if (text=="&Save TOPPAS pipeline")
      {
        bool show = ts && ts->wasChanged();
        actions[i]->setEnabled(show);
      }
      else if (text=="Save TOPPAS pipeline &As")
      {
        bool show = ts && ts->wasChanged();
        actions[i]->setEnabled(show);
      }
      else if (text=="Refresh TOPPAS &parameters") // pipeline menu
      {
        bool show = ts && !(ts->isPipelineRunning());
        actions[i]->setEnabled(show);
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
        if ( (getActive2DWidget()!=0 || getActive3DWidget()!=0) && getActiveCanvas()!=0 )
    		{
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#ffffff;100,#000000");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (*it=="@bg")
    	{
        if ( (getActive2DWidget()!=0 || getActive3DWidget()!=0) && getActiveCanvas()!=0 )
    		{
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#dddddd;100,#000000");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (*it=="@b")
    	{
        if ( (getActive2DWidget()!=0 || getActive3DWidget()!=0) && getActiveCanvas()!=0 )
    		{
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#000000;100,#000000");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (*it=="@r")
    	{
        if ( (getActive2DWidget()!=0 || getActive3DWidget()!=0) && getActiveCanvas()!=0 )
    		{
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#ff0000;100,#ff0000");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (*it=="@g")
    	{
        if ( (getActive2DWidget()!=0 || getActive3DWidget()!=0) && getActiveCanvas()!=0 )
    		{
          Param tmp = getActiveCanvas()->getCurrentLayer().param;
    			tmp.setValue("dot:gradient", "Linear|0,#00ff00;100,#00ff00");
          getActiveCanvas()->setCurrentLayerParameters(tmp);
    		}
    	}
    	else if (*it=="@m")
    	{
        if ( (getActive2DWidget()!=0 || getActive3DWidget()!=0) && getActiveCanvas()!=0 )
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
        addDataFile(*it, false, true);  // add data file but don't show options
    	}
    	else
    	{
    		splash_screen->showMessage((String("Loading file: ") + *it).toQString());
    		splash_screen->repaint();
    		QApplication::processEvents();
    		last_was_plus = false;
        addDataFile(*it, false, true,"",getActiveSpectrumWidget()->getWindowId());
    	}
    }
  }

  void TOPPViewBase::showLogMessage_(TOPPViewBase::LogState state, const String& heading, const String& body)
  {
		//Compose current time string
		DateTime d = DateTime::now();

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

	void TOPPViewBase::showPreferences()
	{
    getActiveCanvas()->showCurrentLayerPreferences();
	}

	void TOPPViewBase::metadataFileDialog()
	{
	 	QStringList files = getFileList_();
		FileHandler fh;
		fh.getOptions().setMetadataOnly(true);
		for (QStringList::iterator it = files.begin(); it!=files.end(); ++it)
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
					dlg.add(exp);
			 	 	dlg.exec();
				}
			}
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

  		//set wait cursor
  		setCursor(Qt::WaitCursor);

			//determine where to copy the data
			UInt new_id = 0;
      if (id != -1) new_id = id;

			if (source == layer_manager_)
			{
				//only the selected row can be dragged => the source layer is the selected layer
        const LayerData& layer = getActiveCanvas()->getCurrentLayer();

         //attach feature, consensus and peak data
        FeatureMapSharedPtrType features = layer.getFeatureMap();
        ExperimentSharedPtrType peaks = layer.getPeakData();
        ConsensusMapSharedPtrType consensus = layer.getConsensusMap();
				vector<PeptideIdentification> peptides = layer.peptides;

				//add the data
        addData(features, consensus, peptides, peaks, layer.type, false, false, true, layer.filename, layer.name, new_id);
			}
      else if (source == spectra_view_treewidget)
			{
        const LayerData& layer = getActiveCanvas()->getCurrentLayer();
        QTreeWidgetItem* item = spectra_view_treewidget->currentItem();
				if (item != 0)
				{
					Size index = (Size)(item->text(3).toInt());
          const ExperimentType::SpectrumType spectrum = (*layer.getPeakData())[index];
					ExperimentType new_exp;
					new_exp.push_back(spectrum);
          ExperimentSharedPtrType new_exp_sptr(new ExperimentType(new_exp));
          FeatureMapSharedPtrType f_dummy(new FeatureMapType());
          ConsensusMapSharedPtrType c_dummy(new ConsensusMapType());
					vector<PeptideIdentification> p_dummy;
          addData(f_dummy, c_dummy, p_dummy, new_exp_sptr, LayerData::DT_CHROMATOGRAM, false, false, true, layer.filename, layer.name, new_id);
				}
			}
			else if (source == 0)
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
    catch(Exception::BaseException& e)
    {
    	showLogMessage_(LS_ERROR,"Error while creating layer",e.what());
    }

		//reset cursor
  	setCursor(Qt::ArrowCursor);
	}


	void TOPPViewBase::updateCurrentPath()
	{
		//do not update if the user disabled this feature.
		if (param_.getValue("preferences:default_path_current")!="true") return;

		//reset
		current_path_ = param_.getValue("preferences:default_path");

		//update if the current layer has a path associated
    if (getActiveCanvas() && getActiveCanvas()->getLayerCount()!=0 && getActiveCanvas()->getCurrentLayer().filename!="")
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

    QWidgetList wl = ws_->windowList();

    // iterate over all windows and determine which need an update
    std::vector<std::pair<const SpectrumWidget *, Size> > needs_update;
    for(int i=0; i!=ws_->windowList().count(); ++i)
    {
      //std::cout << "Number of windows: " << ws_->windowList().count() << std::endl;
      QWidget* w = wl[i];
      const SpectrumWidget* sw = qobject_cast<const SpectrumWidget*>(w);
      if (sw!=0)
      {
        Size lc = sw->canvas()->getLayerCount();

        // determine if widget stores one or more layers for the given filename (->needs update)
        for (Size j=0; j!= lc; ++j)
        {
          //std::cout << "Layer filename: " << sw->canvas()->getLayer(j).filename << std::endl;
          const LayerData& ld = sw->canvas()->getLayer(j);
          if (ld.filename == filename)
          {
            needs_update.push_back(std::pair<const SpectrumWidget *, Size>(sw,j));
          }
        }
      }
    }

    if (needs_update.empty()) // no layer references data of filename
    {
      watcher_->removeFile(filename);  // remove watcher
      return;
    } else if ( !needs_update.empty() )  // at least one layer references data of filename
    {
      //std::cout << "Number of Layers that need update: " << needs_update.size() << std::endl;
      pair<const SpectrumWidget *, Size>& slp = needs_update[0];
      const SpectrumWidget * sw = slp.first;
      Size layer_index = slp.second;

      bool user_wants_update = false;
      if ((String)(param_.getValue("preferences:on_file_change"))=="update automatically") //automatically update
      {
        user_wants_update = true;
      }
      else if ((String)(param_.getValue("preferences:on_file_change"))=="ask") //ask the user if the layer should be updated
      {
        if (watcher_msgbox_==true)
        { // we already have a dialog for that opened... do not ask again
          return;
        }
        // track that we will show the msgbox and we do not need to show it again if file changes once more and the dialog is still open
        watcher_msgbox_=true;
        QMessageBox msg_box;
        QAbstractButton* ok = msg_box.addButton(QMessageBox::Ok);
        msg_box.addButton(QMessageBox::Cancel);
        msg_box.setWindowTitle("Layer data changed");
        msg_box.setText((String("The data of file '") + filename + "' has changed.<BR>Update layers?").toQString());
        msg_box.exec();
        watcher_msgbox_=false;
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
        const LayerData& layer = sw->canvas()->getLayer(layer_index);
        // reload data
        if (layer.type==LayerData::DT_PEAK) //peak data
        {
          try
          {
            FileHandler().loadExperiment(layer.filename,*layer.getPeakData());
          }
          catch(Exception::BaseException& e)
          {
            QMessageBox::critical(this,"Error",(String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
            layer.getPeakData()->clear(true);
          }
          layer.getPeakData()->sortSpectra(true);
          layer.getPeakData()->updateRanges(1);
        }
        else if (layer.type==LayerData::DT_FEATURE) //feature data
        {
          try
          {
            FileHandler().loadFeatures(layer.filename,*layer.getFeatureMap());
          }
          catch(Exception::BaseException& e)
          {
            QMessageBox::critical(this,"Error",(String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
            layer.getFeatureMap()->clear(true);
          }
          layer.getFeatureMap()->updateRanges();
        }
        else if (layer.type==LayerData::DT_CONSENSUS)  //consensus feature data
        {
          try
          {
            ConsensusXMLFile().load(layer.filename,*layer.getConsensusMap());
          }
          catch(Exception::BaseException& e)
          {
            QMessageBox::critical(this,"Error",(String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
            layer.getConsensusMap()->clear(true);
          }
          layer.getConsensusMap()->updateRanges();
        }
        else if (layer.type==LayerData::DT_CHROMATOGRAM) //chromatgram
        {
          //TODO CHROM
          try
          {
            FileHandler().loadExperiment(layer.filename,*layer.getPeakData());
          }
          catch(Exception::BaseException& e)
          {
            QMessageBox::critical(this,"Error",(String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
            layer.getPeakData()->clear(true);
          }
          layer.getPeakData()->sortChromatograms(true);
          layer.getPeakData()->updateRanges(1);

        }
        /*      else if (layer.type == LayerData::DT_IDENT) // identifications
      {
        try
        {
          vector<ProteinIdentification> proteins;
          IdXMLFile().load(layer.filename, proteins, layer.peptides);
        }
        catch(Exception::BaseException& e)
        {
          QMessageBox::critical(this,"Error",(String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
          layer.peptides.clear();
        }
      }
*/
      }

      // update all layers that need an update
      for (Size i=0; i!= needs_update.size(); ++i)
      {
        pair<const SpectrumWidget *, Size>& slp = needs_update[i];
        const SpectrumWidget * sw = slp.first;
        Size layer_index = slp.second;
        sw->canvas()->updateLayer(layer_index);
      }
    }
    /*
    {
      //update the layer if the user choosed to do so
      if (update)
      {
        emit sendStatusMessage(String("Updating layer '") + getLayer(j).name + "' (file changed).",0);
        updateLayer_(j);
        emit sendStatusMessage(String("Finished updating layer '") + getLayer(j).name + "'.",5000);
      }
    }
    */
    updateLayerBar();
    updateViewBar();
    updateFilterBar();
    updateMenu();

    // temporarily remove and read filename from watcher_ as a workaround for bug #233
    watcher_->removeFile(filename);
    watcher_->addFile(filename);
  }

  void TOPPViewBase::setTOPPASTabEnabled(bool enabled)
  {
    if (views_tabwidget_->isTabEnabled(2) != enabled)
    {
      views_tabwidget_->setTabEnabled(2, enabled);
    }
  }

  // ******************************************TOPPAS events**********************************
  void TOPPViewBase::newPipeline()
  {
    TOPPASWidget* tw = new TOPPASWidget(Param(), ws_, toppas_tmp_path_);
    TOPPASScene* ts = tw->getScene();
    connect (ts, SIGNAL(selectionCopied(TOPPASScene*)), this, SLOT(saveToClipboard(TOPPASScene*)));
    connect (ts, SIGNAL(requestClipboardContent()), this, SLOT(sendClipboardContent()));
    connect (ts, SIGNAL(saveMe()), this, SLOT(savePipeline()));
    connect (ts, SIGNAL(mainWindowNeedsUpdate()), this, SLOT(updateMenu()));
    showTOPPipelineInWindow_(tw, "(Untitled)");
  }

  void TOPPViewBase::includePipeline()
  {
    QString file_name = QFileDialog::getOpenFileName(this, tr("Include workflow"), current_path_.toQString(), tr("TOPPAS pipelines (*.toppas)"));
    addTOPPASFile(file_name, false);
  }

  void TOPPViewBase::savePipeline()
  {
    TOPPASWidget* w = 0;
    QObject* sendr = QObject::sender();
    QAction* save_button_clicked = qobject_cast<QAction*>(sendr);

    if (!save_button_clicked)
    {
      // scene has requested to be saved
      TOPPASScene* ts = qobject_cast<TOPPASScene*>(sendr);
      if (ts && ts->views().size() > 0)
      {
        w = qobject_cast<TOPPASWidget*>(ts->views().first());
      }
    }
    else
    {
      w = getActiveTOPPASWidget();
    }

    if (!w)
    {
      return;
    }

    String file_name = w->getScene()->getSaveFileName();
    if (file_name != "")
    {
      if (!file_name.hasSuffix(".toppas"))
      {
        file_name += ".toppas";
      }
      w->getScene()->store(file_name);
    }
    else
    {
      TOPPASBase::savePipelineAs(w, current_path_.toQString());
    }
  }

  void TOPPViewBase::saveCurrentPipelineAs()
  {
    TOPPASWidget* w = getActiveTOPPASWidget();
    QString file_name = TOPPASBase::savePipelineAs(w, current_path_.toQString());
    if (file_name != "")
    {
      QString caption = File::basename(file_name).toQString();
      tab_bar_->setTabText(tab_bar_->currentIndex(), caption);
    }
  }

  void TOPPViewBase::loadPipelineResourceFile()
  {
    TOPPASWidget* w = getActiveTOPPASWidget();
    TOPPASBase::loadPipelineResourceFile(w, current_path_.toQString());
  }

  void TOPPViewBase::savePipelineResourceFile()
  {
    TOPPASWidget* w = getActiveTOPPASWidget();
    TOPPASBase::savePipelineResourceFile(w, current_path_.toQString());
  }

  void TOPPViewBase::refreshPipelineParameters()
  {
    TOPPASWidget* w = getActiveTOPPASWidget();
    QString file_name = TOPPASBase::refreshPipelineParameters(w, current_path_.toQString());
    if (file_name != "")
    {
      QString caption = File::basename(file_name).toQString();
      tab_bar_->setTabText(tab_bar_->currentIndex(), caption);
    }
  }

  void TOPPViewBase::runPipeline()
  {
    TOPPASWidget* w = getActiveTOPPASWidget();
    if (w)
    {
      w->getScene()->runPipeline();
    }
  }

  void TOPPViewBase::abortPipeline()
  {
    TOPPASWidget* w = getActiveTOPPASWidget();
    if (w)
    {
      w->getScene()->abortPipeline();
    }
    updateMenu();
  }

  void TOPPViewBase::showPipelineFinishedLogMessage()
  {
    showLogMessage_(LS_NOTICE, "Entire pipeline execution finished!", "");
  }

  void TOPPViewBase::saveToClipboard(TOPPASScene* scene)
  {
    if (toppas_clipboard_scene_ != 0)
    {
      delete toppas_clipboard_scene_;
      toppas_clipboard_scene_ = 0;
    }
    toppas_clipboard_scene_ = scene;
  }

  void TOPPViewBase::sendClipboardContent()
  {
    TOPPASScene* sndr = qobject_cast<TOPPASScene*>(QObject::sender());
    if (sndr != 0)
    {
      sndr->setClipboard(toppas_clipboard_scene_);
    }
  }

  void TOPPViewBase::toolStarted()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (type != "")
      {
        text += " ("+type+")";
      }
      text += " started. Processing ...";

      showLogMessage_(LS_NOTICE, text, "");
    }
    updateMenu();
  }

  void TOPPViewBase::toolFinished()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (type != "")
      {
        text += " ("+type+")";
      }
      text += " finished!";

      showLogMessage_(LS_NOTICE, text, "");
    }
    updateMenu();
  }

  void TOPPViewBase::toolCrashed()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (type != "")
      {
        text += " ("+type+")";
      }
      text += " crashed!";

      showLogMessage_(LS_ERROR, text, "");
    }
    updateMenu();
  }

  void TOPPViewBase::toolFailed()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (type != "")
      {
        text += " ("+type+")";
      }
      text += " failed!";

      showLogMessage_(LS_ERROR, text, "");
    }
    updateMenu();
  }

  void TOPPViewBase::outputVertexFinished(const String& file)
  {
    String text = "Output file '" + file + "' written.";
    showLogMessage_(LS_NOTICE, text, "");
  }

  void TOPPViewBase::updateTOPPOutputLog(const QString& out)
  {
    QString text = out; // shortened version for now (if we reintroduce simultaneous tool execution,
                        // we need to rethink this (probably only trigger this slot when tool 100% finished)

    //show log if there is output
    qobject_cast<QWidget*>(log_->parent())->show();

    //update log_
    log_->moveCursor(QTextCursor::End,QTextCursor::MoveAnchor); // move cursor to end, since text is inserted at cursor
		log_->insertPlainText(text);
  }

  void TOPPViewBase::openFilesInTOPPView(QStringList files)
  {
    foreach(QString s, files)
    {
      addDataFile(s, false, false, s);
    }
  }

  //

  TOPPViewBase::~TOPPViewBase()
  {
    savePreferences();
    abortTOPPTool();

    // dispose behavior
    if(identificationview_behavior_ != 0)
    {
      delete(identificationview_behavior_);
    }

    if(spectraview_behavior_ != 0)
    {
      delete(spectraview_behavior_);
    }
  }

} //namespace OpenMS
