// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPViewMenu.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/RecentFilesMenu.h>

#include <QAction>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>

using namespace std;

namespace OpenMS
{

  FS_TV operator+(const TV_STATUS left, const TV_STATUS right)
  {
    FS_TV r(left);
    r += right;
    return r;
  }

  FS_LAYER OPENMS_GUI_DLLAPI operator+(const LayerDataBase::DataType left, const LayerDataBase::DataType right)
  {
    FS_LAYER r;
    r += left;
    r += right;
    return r;
  }

  TOPPViewMenu::TOPPViewMenu(TOPPViewBase* const parent, EnhancedWorkspace* const ws, RecentFilesMenu* const recent_files)
    : QObject()
      //parent_(parent)
  {
    QAction* action; ///< for adding tool tips to actions
    
    QMenu* m_file = new QMenu("&File", parent);
    m_file->setToolTipsVisible(true);
    parent->menuBar()->addMenu(m_file);

    // we explicitly pass an empty Path here using a Lambda, since using the default `..., parent, &TOPPViewBase::openFilesByDialog, ...`)
    // passes a "0" as argument (Qt bug?)
    m_file->addAction("&Open file", parent, [parent]() { parent->openFilesByDialog(""); })->setShortcut(Qt::CTRL | Qt::Key_O);
    m_file->addAction("Open &example file", parent, [parent]() { parent->openFilesByDialog(File::getOpenMSDataPath() + "/examples/"); })->setShortcut(Qt::CTRL | Qt::Key_E);
    action = m_file->addAction("&Close tab", parent, &TOPPViewBase::closeTab);
    action->setShortcut(Qt::CTRL | Qt::Key_W);
    addAction_(action, TV_STATUS::HAS_CANVAS);
    m_file->addSeparator();

    // Meta data
    action = m_file->addAction("&Show meta data (file)", parent, &TOPPViewBase::metadataFileDialog);
    action->setToolTip("Load a file's meta information without actually loading the data.");

    m_file->addSeparator();

    // Recent files
    m_file->addMenu(recent_files->getMenu()); // updates automatically via RecentFilesMenu class, since this is just a pointer

    m_file->addSeparator();
    
    // Specifically set the role of the Preferences item. Additionally we have to avoid adding other action items that are
    // called preferences/config/options and have the default TextHeuristicRole because otherwise they will overwrite the macOS specific
    // menu entry under Application -> Preferences...
    // m_file->addAction("&Preferences", parent, &TOPPViewBase::preferencesDialog);
    auto pref = new QAction("&Preferences", parent);
    pref->setMenuRole(QAction::PreferencesRole);
    pref->setEnabled(true);
    m_file->addAction(pref);
    connect(pref, &QAction::triggered, parent, &TOPPViewBase::preferencesDialog);
      
    m_file->addAction("&Quit", qApp, SLOT(quit()));

    // Tools menu
    QMenu* m_tools = new QMenu("&Tools", parent);
    m_tools->setToolTipsVisible(true);
    parent->menuBar()->addMenu(m_tools);
    action = m_tools->addAction("&Select data range", parent, &TOPPViewBase::showGoToDialog);
    action->setShortcut(Qt::CTRL | Qt::Key_G);
    addAction_(action, TV_STATUS::HAS_LAYER);
    action = m_tools->addAction("&Edit meta data", parent, &TOPPViewBase::editMetadata);
    action->setShortcut(Qt::CTRL | Qt::Key_M);
    addAction_(action, TV_STATUS::HAS_LAYER);
    addAction_(m_tools->addAction("&Statistics", parent, &TOPPViewBase::layerStatistics),
      TV_STATUS::HAS_LAYER);
    m_tools->addSeparator();
    action = m_tools->addAction("Apply TOPP tool (whole layer)", parent, &TOPPViewBase::showTOPPDialog);
    action->setShortcut(Qt::CTRL | Qt::Key_T);
    action->setData(false);
    addAction_(action, TV_STATUS::HAS_LAYER + TV_STATUS::TOPP_IDLE);
    action = m_tools->addAction("Apply TOPP tool (visible layer data)", parent, &TOPPViewBase::showTOPPDialog);
    action->setShortcut(Qt::CTRL | Qt::SHIFT | Qt::Key_T);
    action->setData(true);
    addAction_(action, TV_STATUS::HAS_LAYER + TV_STATUS::TOPP_IDLE);
    action = m_tools->addAction("Rerun TOPP tool", parent, &TOPPViewBase::rerunTOPPTool);
    action->setShortcut(Qt::Key_F4);
    addAction_(action, TV_STATUS::HAS_LAYER + TV_STATUS::TOPP_IDLE);
    m_tools->addSeparator();

    action = m_tools->addAction("&Annotate with AccurateMassSearch results", parent, &TOPPViewBase::annotateWithAMS);
    action->setShortcut(Qt::CTRL | Qt::Key_A);
    action->setToolTip("Annotate Peak layer with a featureXML from the AccurateMassSearch tool");
    addAction_(action, TV_STATUS::HAS_LAYER, FS_LAYER(LayerDataBase::DT_PEAK));

    action = m_tools->addAction("&Annotate with peptide identifications", parent, &TOPPViewBase::annotateWithID);
    action->setShortcut(Qt::CTRL | Qt::Key_I);
    action->setToolTip("Annotate a Peak or Feature or Consensus layer with peptide identifications");
    addAction_(action, TV_STATUS::HAS_LAYER, LayerDataBase::DT_PEAK + LayerDataBase::DT_FEATURE + LayerDataBase::DT_CONSENSUS);

    action = m_tools->addAction("&Annotate with OpenSwath transitions", parent, &TOPPViewBase::annotateWithOSW);
    action->setShortcut(Qt::CTRL | Qt::Key_P);
    action->setToolTip("Annotate Chromatogram layer with OSW transition id data from OpenSwathWorkflow or pyProphet");
    addAction_(action, TV_STATUS::HAS_LAYER, FS_LAYER(LayerDataBase::DT_CHROMATOGRAM));
    
    action = addAction_(m_tools->addAction("Align spectra", parent, &TOPPViewBase::showSpectrumAlignmentDialog),
      TV_STATUS::HAS_MIRROR_MODE);
    action->setToolTip("Only available in 1D View for mirrored (flipped) spectra. To flip, use the Layer View and right click a layer.");
    
    m_tools->addAction("Generate theoretical spectrum", parent, &TOPPViewBase::showSpectrumGenerationDialog);

    // Layer menu
    QMenu* m_layer = new QMenu("&Layer", parent);
    m_layer->setToolTipsVisible(true);
    parent->menuBar()->addMenu(m_layer);
    action = m_layer->addAction("Save all data", parent, &TOPPViewBase::saveLayerAll);
    action->setShortcut(Qt::CTRL | Qt::Key_S);
    addAction_(action, TV_STATUS::HAS_LAYER);
    action = m_layer->addAction("Save visible data", parent, &TOPPViewBase::saveLayerVisible);
    action->setShortcut(Qt::CTRL | Qt::SHIFT | Qt::Key_S);
    addAction_(action, TV_STATUS::HAS_LAYER);
    m_layer->addSeparator();
    action = m_layer->addAction("Show/hide grid lines", parent, &TOPPViewBase::toggleGridLines);
    action->setShortcut(Qt::CTRL | Qt::Key_R);
    addAction_(action, TV_STATUS::HAS_LAYER);
    action = m_layer->addAction("Show/hide axis legends", parent, &TOPPViewBase::toggleAxisLegends);
    action->setShortcut(Qt::CTRL | Qt::Key_L);
    addAction_(action, TV_STATUS::HAS_CANVAS);
    action = addAction_(m_layer->addAction("Show/hide automated m/z annotations", parent, &TOPPViewBase::toggleInterestingMZs),
      TV_STATUS::IS_1D_VIEW);
    action->setToolTip("Only available in 1D View");
    m_layer->addSeparator();
    
    // Do not call it preferences without disabling text heuristics role.
    addAction_(m_layer->addAction("Layer preferences", parent, &TOPPViewBase::showPreferences),
      TV_STATUS::HAS_LAYER);

    // Windows menu
    m_windows_ = new QMenu("&Windows", parent);
    m_windows_->setToolTipsVisible(true);
    parent->menuBar()->addMenu(m_windows_);
    m_windows_->addAction("&Cascade", ws, &EnhancedWorkspace::cascadeSubWindows);
    m_windows_->addAction("&Tile automatic", ws, &EnhancedWorkspace::tileSubWindows);
    m_windows_->addAction(QIcon(":/tile_vertical.png"), "Tile &vertical", ws, &EnhancedWorkspace::tileVertical);
    m_windows_->addAction(QIcon(":/tile_horizontal.png"), "Tile &horizontal", ws, &EnhancedWorkspace::tileHorizontal);
    // link / unlink
    action = m_windows_->addAction("Link/Unlink &Zoom", parent, &TOPPViewBase::linkZoom);
    action->setToolTip("Zoom all open tab windows to the same coordinates concurrently (requires the same view dimension; e.g. all 2D views will show the same RT/mz windows). Most effective when used in tiled Windows view (see Windows -> tiling)");
    m_windows_->addSeparator();

    // Help menu
    QMenu* m_help = new QMenu("&Help", parent);
    m_help->setToolTipsVisible(true);
    parent->menuBar()->addMenu(m_help);
    m_help->addAction(QWhatsThis::createAction(m_help));
    m_help->addSeparator();
    m_help->addAction("OpenMS website", []() { GUIHelpers::openURL("http://www.OpenMS.de"); });
    m_help->addAction("Tutorials and documentation", []() { GUIHelpers::openURL("html/index.html"); })->setShortcut(Qt::Key_F1);

    m_help->addSeparator();

    // Note: it is important to pass parent by value, since the lambda will be evaluated later, 
    // even after this function returned and parent reference would be out of scope.
    m_help->addAction("&About", [parent]() {QApplicationTOPP::showAboutDialog(parent, "TOPPView"); });
  }

  void TOPPViewMenu::update(const FS_TV status, const LayerDataBase::DataType layer_type)
  {
    for (auto& ar : menu_items_)
    { // only disable if not supported by the view. This way, the user can still see the item (greyed out) and its ToolTip (for how to activate the item)
      ar.enableAction(status, layer_type);
    }
  }

  void TOPPViewMenu::addWindowToggle(QAction* const window_toggle)
  {
    m_windows_->addAction(window_toggle);
  }

  QAction* TOPPViewMenu::addAction_(QAction* action, const TV_STATUS req, const FS_LAYER type)
  {
    menu_items_.emplace_back(action, req, type);
    return action;
  }
  QAction* TOPPViewMenu::addAction_(QAction* action, const FS_TV req, const FS_LAYER type)
  {
    menu_items_.emplace_back(action, req, type);
    return action;
  }


  void TOPPViewMenu::ActionRequirement_::enableAction(const FS_TV status, const LayerDataBase::DataType layer_type)
  {
    bool status_ok = status.isSuperSetOf(needs_);
    bool layer_ok = layer_set_.isSuperSetOf(layer_type) || layer_set_.empty();
    this->action_->setEnabled(status_ok && layer_ok);
  }

} //Namespace
