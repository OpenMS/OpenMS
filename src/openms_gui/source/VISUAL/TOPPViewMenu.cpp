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

  TOPPViewMenu::TOPPViewMenu(TOPPViewBase* const parent, EnhancedWorkspace* const ws, RecentFilesMenu* const recent_files)
    : QObject()
      //parent_(parent)
  {
    QAction* action; ///< for adding tool tips to actions
    
    QMenu* m_file = new QMenu("&File", parent);
    m_file->setToolTipsVisible(true);
    parent->menuBar()->addMenu(m_file);

    m_file->addAction("&Open file", parent, &TOPPViewBase::openFilesByDialog, Qt::CTRL + Qt::Key_O);
    m_file->addAction("Open &example file", parent, [parent]() { parent->openFilesByDialog(File::getOpenMSDataPath() + "/examples/"); }, Qt::CTRL + Qt::Key_E);
    addAction_(m_file->addAction("&Close tab", parent, &TOPPViewBase::closeTab, Qt::CTRL + Qt::Key_W),
               TV_STATUS::HAS_CANVAS);
    m_file->addSeparator();

    // Meta data
    action = m_file->addAction("&Show meta data (file)", parent, &TOPPViewBase::metadataFileDialog);
    action->setToolTip("Load a file's meta information without actually loading the data.");

    m_file->addSeparator();

    // Recent files
    m_file->addMenu(recent_files->getMenu()); // updates automatically via RecentFilesMenu class, since this is just a pointer

    m_file->addSeparator();
    addAction_(m_file->addAction("&Preferences", parent, &TOPPViewBase::preferencesDialog),
               TV_STATUS::HAS_LAYER);
    m_file->addAction("&Quit", qApp, SLOT(quit()));

    // Tools menu
    QMenu* m_tools = new QMenu("&Tools", parent);
    m_tools->setToolTipsVisible(true);
    parent->menuBar()->addMenu(m_tools);
    addAction_(m_tools->addAction("&Select data range", parent, &TOPPViewBase::showGoToDialog, Qt::CTRL + Qt::Key_G),
      TV_STATUS::HAS_LAYER);
    addAction_(m_tools->addAction("&Edit meta data", parent, &TOPPViewBase::editMetadata, Qt::CTRL + Qt::Key_M),
      TV_STATUS::HAS_LAYER);
    addAction_(m_tools->addAction("&Statistics", parent, &TOPPViewBase::layerStatistics),
      TV_STATUS::HAS_LAYER);
    m_tools->addSeparator();
    auto item = addAction_(m_tools->addAction("Apply TOPP tool (whole layer)", parent, &TOPPViewBase::showTOPPDialog, Qt::CTRL + Qt::Key_T),
        TV_STATUS::HAS_LAYER + TV_STATUS::TOPP_IDLE);
    item.action->setData(false);
    item = addAction_(m_tools->addAction("Apply TOPP tool (visible layer data)", parent, &TOPPViewBase::showTOPPDialog, Qt::CTRL + Qt::SHIFT + Qt::Key_T),
      TV_STATUS::HAS_LAYER + TV_STATUS::TOPP_IDLE);
    item.action->setData(true);
    addAction_(m_tools->addAction("Rerun TOPP tool", parent, &TOPPViewBase::rerunTOPPTool, Qt::Key_F4),
      TV_STATUS::HAS_LAYER + TV_STATUS::TOPP_IDLE);
    m_tools->addSeparator();
    addAction_(m_tools->addAction("&Annotate with identification", parent, &TOPPViewBase::annotateWithID, Qt::CTRL + Qt::Key_I),
      TV_STATUS::HAS_LAYER);
    action = addAction_(m_tools->addAction("Align spectra", parent, &TOPPViewBase::showSpectrumAlignmentDialog),
      TV_STATUS::HAS_MIRROR_MODE).action;
    action->setToolTip("Only available in 1D View for mirrored (flipped) spectra. To flip, use the Layer View and right click a layer.");
    m_tools->addAction("Generate theoretical spectrum", parent, &TOPPViewBase::showSpectrumGenerationDialog);

    // Layer menu
    QMenu* m_layer = new QMenu("&Layer", parent);
    m_layer->setToolTipsVisible(true);
    parent->menuBar()->addMenu(m_layer);
    addAction_(m_layer->addAction("Save all data", parent, &TOPPViewBase::saveLayerAll, Qt::CTRL + Qt::Key_S),
      TV_STATUS::HAS_LAYER);
    addAction_(m_layer->addAction("Save visible data", parent, &TOPPViewBase::saveLayerVisible, Qt::CTRL + Qt::SHIFT + Qt::Key_S),
      TV_STATUS::HAS_LAYER);
    m_layer->addSeparator();
    addAction_(m_layer->addAction("Show/hide grid lines", parent, &TOPPViewBase::toggleGridLines, Qt::CTRL + Qt::Key_R),
      TV_STATUS::HAS_LAYER);
    addAction_(m_layer->addAction("Show/hide axis legends", parent, &TOPPViewBase::toggleAxisLegends, Qt::CTRL + Qt::Key_L),
      TV_STATUS::HAS_CANVAS);
    action = addAction_(m_layer->addAction("Show/hide automated m/z annotations", parent, &TOPPViewBase::toggleInterestingMZs),
      TV_STATUS::IS_1D_VIEW).action;
    action->setToolTip("Only available in 1D View");
    m_layer->addSeparator();
    addAction_(m_layer->addAction("Preferences", parent, &TOPPViewBase::showPreferences),
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
    m_help->addAction("OpenMS website", [&]() { GUIHelpers::openURL("http://www.OpenMS.de"); });
    m_help->addAction("Tutorials and documentation", [&]() { GUIHelpers::openURL("html/index.html"); }, Qt::Key_F1);

    m_help->addSeparator();
    m_help->addAction("&About", [&]() {QApplicationTOPP::showAboutDialog(parent, "TOPPView"); });
  }

  void TOPPViewMenu::update(const FS_TV status)
  {
    for (auto& ar : menu_items_)
    { // only disable if not supported by the view. This way, the user can still see the item (greyed out) and its ToolTip (for how to activate the item)
      ar.action->setEnabled(status.isSuperSetOf(ar.needs));
    }
  }

  void TOPPViewMenu::addWindowToggle(QAction* const window_toggle)
  {
    m_windows_->addAction(window_toggle);
  }

  const TOPPViewMenu::ActionRequirement_& TOPPViewMenu::addAction_(QAction* action, const TV_STATUS req)
  {
    menu_items_.emplace_back(action, req);
    return menu_items_.back();
  }
  const TOPPViewMenu::ActionRequirement_& TOPPViewMenu::addAction_(QAction* action, const FS_TV req)
  {
    menu_items_.emplace_back(action, req);
    return menu_items_.back();
  }


} //Namespace
