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

#include <OpenMS/VISUAL/LayerListView.h>

#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum3DWidget.h>

#include <QtWidgets/QListWidgetItem>

using namespace std;

namespace OpenMS
{
  LayerListView::LayerListView(QWidget* parent)
    : QListWidget(parent)
  {
    const auto help = "Layer bar<BR>"
                       "<BR>Here the available layers are shown. Left-click on a layer to select it."
                       "<BR>Layers can be shown and hidden using the checkboxes in front of the name."
                       "<BR>Renaming and removing a layer is possible through the context menu."
                       "<BR>Dragging a layer to the tab bar copies the layer."
                       "<BR>Double-clicking a layer open its preferences."
                       "<BR>You can use the 'PageUp' and 'PageDown' buttons to change the selected layer.";
    this->setWhatsThis(help);
    this->setToolTip(help);

    setDragEnabled(true);
    connect(this, &LayerListView::currentRowChanged, this, &LayerListView::currentRowChangedAction_);
    connect(this, &LayerListView::itemChanged, this, &LayerListView::itemChangedAction_);
    connect(this, &QListWidget::itemDoubleClicked, this, &LayerListView::itemDoubleClickedAction_);
  }

  void LayerListView::update(SpectrumWidget* active_widget)
  {
    // reset items
    this->clear();

    spectrum_widget_ = active_widget;
    // during program exit, this could be called after SpectrumWidgets are gone
    if (spectrum_widget_ == nullptr) return;

    SpectrumCanvas* cc = spectrum_widget_->canvas();
    if (cc == nullptr) return;

    // determine if this is a 1D view (for text color)
    bool is_1d_view = (dynamic_cast<Spectrum1DCanvas*>(cc) != nullptr);

    this->blockSignals(true);
    RAIICleanup cl([&]() { this->blockSignals(false); });

    for (Size i = 0; i < cc->getLayerCount(); ++i)
    {
      const LayerData& layer = cc->getLayer(i);

      // add item
      QListWidgetItem* item = new QListWidgetItem(this);
      QString name = layer.getDecoratedName().toQString();

      item->setText(name);
      item->setToolTip(layer.filename.toQString());

      if (is_1d_view)
      {
        QPixmap icon(7, 7);
        icon.fill(QColor(layer.param.getValue("peak_color").toQString()));
        item->setIcon(icon);
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
      if (i == cc->getCurrentLayerIndex())
      {
        this->setCurrentItem(item);
      }
    }
  }

  void LayerListView::currentRowChangedAction_(int i)
  {
    // after adding a layer, i is -1. TODO: check if this is the correct behaviour
    if (i != -1)
    {
      spectrum_widget_->canvas()->activateLayer(i); // emits layerActivated in TOPPView
    }
  }
  void LayerListView::itemChangedAction_(QListWidgetItem* item)
  {
    int layer = this->row(item);
    bool visible = spectrum_widget_->canvas()->getLayer(layer).visible;

    if (item->checkState() == Qt::Unchecked && visible)
    {
      spectrum_widget_->canvas()->changeVisibility(layer, false);
      emit layerDataChanged();
    }
    else if (item->checkState() == Qt::Checked && !visible)
    {
      spectrum_widget_->canvas()->changeVisibility(layer, true);
      emit layerDataChanged();
    }
  }
  void LayerListView::contextMenuEvent(QContextMenuEvent* event)
  {
    QListWidgetItem* item = this->itemAt(event->pos());
    if (!item) return;

    int layer_idx = this->row(item);
    QMenu* context_menu = new QMenu(this);
    
    context_menu->addAction("Rename", [&]() {
      QString name = QInputDialog::getText(this, "Rename layer", "Name:", QLineEdit::Normal, spectrum_widget_->canvas()->getLayerName(layer_idx).toQString());
      if (name != "")
      {
        spectrum_widget_->canvas()->setLayerName(layer_idx, name);
        emit layerDataChanged();
      }});

    context_menu->addAction("Delete", [&]() {
      spectrum_widget_->canvas()->removeLayer(layer_idx);
      emit layerDataChanged();
    });

    auto widget1D = qobject_cast<Spectrum1DWidget*>(spectrum_widget_);
    if (widget1D != nullptr)
    {
      if (widget1D->canvas()->getLayer(layer_idx).flipped)
      {
        context_menu->addAction("Flip upwards (1D)", [&]() {
          widget1D->canvas()->flipLayer(layer_idx);
          widget1D->canvas()->setMirrorModeActive(widget1D->canvas()->flippedLayersExist());
        });
        emit layerDataChanged();
      }
      else
      {
        context_menu->addAction("Flip downwards (1D)", [&]() {
          widget1D->canvas()->flipLayer(layer_idx);
          widget1D->canvas()->setMirrorModeActive(true);
        });
        emit layerDataChanged();
      }
    }

    context_menu->addSeparator();
    context_menu->addAction("Preferences", [&]() {
      spectrum_widget_->canvas()->showCurrentLayerPreferences();
    });

    context_menu->exec(this->mapToGlobal(event->pos()));
  }

  void LayerListView::itemDoubleClickedAction_(QListWidgetItem* /*item*/)
  {
    spectrum_widget_->canvas()->showCurrentLayerPreferences();
  }
} //Namespace
