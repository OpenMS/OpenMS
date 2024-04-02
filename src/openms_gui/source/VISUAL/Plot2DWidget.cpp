// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/Plot1DWidget.h>
#include <OpenMS/VISUAL/Plot2DWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Plot2DGoToDialog.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>

#include <OpenMS/VISUAL/LayerDataConsensus.h>
#include <OpenMS/VISUAL/LayerDataFeature.h>

#include <QtWidgets/QPushButton>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QCheckBox>
#include <QtCore/QTimer>

using namespace std;

namespace OpenMS
{
  using namespace Internal;
  using namespace Math;

  Plot2DWidget::Plot2DWidget(const Param& preferences, QWidget* parent) :
    PlotWidget(preferences, parent)
  {
    setCanvas_(new Plot2DCanvas(preferences, this), 1, 2);

    // add axes
    y_axis_->setMinimumWidth(50);

    // add projections
    grid_->setColumnStretch(2, 3);
    grid_->setRowStretch(1, 3);

    projection_onto_X_ = new Plot1DWidget(Param(), DIM::Y, this);
    projection_onto_X_->hide();
    grid_->addWidget(projection_onto_X_, 0, 1, 1, 2);

    projection_onto_Y_ = new Plot1DWidget(Param(), DIM::X, this);
    projection_onto_Y_->hide();
    grid_->addWidget(projection_onto_Y_, 1, 3, 2, 1);

    connect(canvas(), &Plot2DCanvas::showProjections, this, &Plot2DWidget::showProjections_);
    connect(canvas(), &Plot2DCanvas::toggleProjections, this, &Plot2DWidget::toggleProjections);
    connect(canvas(), &Plot2DCanvas::visibleAreaChanged, this, &Plot2DWidget::autoUpdateProjections_);
    // delegate signals from canvas
    connect(canvas(), &Plot2DCanvas::showSpectrumAsNew1D, this, &Plot2DWidget::showSpectrumAsNew1D);
    connect(canvas(), &Plot2DCanvas::showChromatogramsAsNew1D, this, &Plot2DWidget::showChromatogramsAsNew1D);
    connect(canvas(), &Plot2DCanvas::showCurrentPeaksAsIonMobility, this, &Plot2DWidget::showCurrentPeaksAsIonMobility);
    connect(canvas(), &Plot2DCanvas::showCurrentPeaksAs3D, this, &Plot2DWidget::showCurrentPeaksAs3D);
    // add projections box
    projection_box_ = new QGroupBox("Projections", this);
    projection_box_->hide();
    grid_->addWidget(projection_box_, 0, 3);
    QGridLayout* box_grid = new QGridLayout(projection_box_);

    QLabel* label = new QLabel("Peaks: ");
    box_grid->addWidget(label, 0, 0);
    projection_peaks_ = new QLabel("");
    box_grid->addWidget(projection_peaks_, 0, 1);

    label = new QLabel("Intensity sum: ");
    box_grid->addWidget(label, 1, 0);
    projection_sum_ = new QLabel("");
    box_grid->addWidget(projection_sum_, 1, 1);

    label = new QLabel("Maximum intensity: ");
    box_grid->addWidget(label, 2, 0);
    projection_max_ = new QLabel("");
    box_grid->addWidget(projection_max_, 2, 1);

    box_grid->setRowStretch(3, 2);

    QPushButton* button = new QPushButton("Update", projection_box_);
    connect(button, &QPushButton::clicked, canvas(), &Plot2DCanvas::pickProjectionLayer);
    box_grid->addWidget(button, 4, 0);

    projections_auto_ = new QCheckBox("Auto-update", projection_box_);
    projections_auto_->setWhatsThis("When activated, projections are automatically updated one second after the last change of the visible area.");
    projections_auto_->setChecked(true);
    box_grid->addWidget(projections_auto_, 4, 1);

    //set up projections auto-update
    projections_timer_ = new QTimer(this);
    projections_timer_->setSingleShot(true);
    projections_timer_->setInterval(1000);
    connect(projections_timer_, &QTimer::timeout, canvas(), &Plot2DCanvas::pickProjectionLayer);
  }

  void Plot2DWidget::projectionInfo_(int peaks, double intensity, double max)
  {
    projection_peaks_->setText(QString::number(peaks));
    projection_sum_->setText(QString::number(intensity, 'f', 1));
    projection_max_->setText(QString::number(max, 'f', 1));
  }

  void Plot2DWidget::recalculateAxes_()
  {
    // set names
    x_axis_->setLegend(string(canvas()->getMapper().getDim(DIM::X).getDimName()));
    y_axis_->setLegend(string(canvas()->getMapper().getDim(DIM::Y).getDimName()));

    const auto& area = canvas()->getVisibleArea().getAreaXY();
    x_axis_->setAxisBounds(area.minX(), area.maxX());
    y_axis_->setAxisBounds(area.minY(), area.maxY());
  }

  void Plot2DWidget::toggleProjections()
  {
    if (projectionsVisible())
    {
      setMinimumSize(250, 250);
      projection_box_->hide();
      projection_onto_Y_->hide();
      projection_onto_X_->hide();
      grid_->setColumnStretch(3, 0);
      grid_->setRowStretch(0, 0);
    }
    else
    {
      setMinimumSize(500, 500);
      canvas()->pickProjectionLayer();
    }
  }

  //  projections
  void Plot2DWidget::showProjections_(const LayerDataBase* source_layer)
  {
    auto [projection_ontoX, projection_ontoY, stats] =
      source_layer->getProjection(canvas_->getMapper().getDim(DIM::X).getUnit(), canvas_->getMapper().getDim(DIM::Y).getUnit(), canvas_->getVisibleArea().getAreaUnit());

    projectionInfo_(stats.number_of_datapoints, stats.sum_intensity, stats.max_intensity);

    auto va = canvas()->getVisibleArea().getAreaXY();

    projection_onto_Y_->showLegend(false);
    projection_onto_Y_->setMapper({{DIM_UNIT::INT, canvas_->getMapper().getDim(DIM::Y).getUnit()}}); // must be done before addLayer()
    projection_onto_Y_->canvas()->removeLayers();
    projection_onto_Y_->canvas()->addLayer(std::move(projection_ontoY));
    // manually set projected unit, since 'addPeakLayer' will guess a visible area, but we want the exact same scaling
    projection_onto_Y_->canvas()->setVisibleAreaY(va.minY(), va.maxY());
    grid_->setColumnStretch(3, 2);

    projection_onto_X_->showLegend(false);
    projection_onto_X_->setMapper({{canvas_->getMapper().getDim(DIM::X).getUnit(), DIM_UNIT::INT}}); // must be done before addLayer()
    projection_onto_X_->canvas()->removeLayers();
    projection_onto_X_->canvas()->addLayer(std::move(projection_ontoX));
    // manually set projected unit, since 'addPeakLayer' will guess a visible area, but we want the exact same scaling
    projection_onto_X_->canvas()->setVisibleAreaX(va.minX(), va.maxX());
    grid_->setRowStretch(0, 2);
    
    projection_box_->show();
    projection_onto_X_->show();
    projection_onto_Y_->show();
  }

  const Plot1DWidget* Plot2DWidget::getProjectionOntoX() const
  {
    return projection_onto_X_;
  }

  const Plot1DWidget* Plot2DWidget::getProjectionOntoY() const
  {
    return projection_onto_Y_;
  }

  void Plot2DWidget::showGoToDialog()
  {
    Plot2DGoToDialog goto_dialog(this,
      canvas_->getMapper().getDim(DIM::X).getDimNameShort(),
      canvas_->getMapper().getDim(DIM::Y).getDimNameShort());
    //set range
    const auto& area = canvas()->getVisibleArea().getAreaXY();
    goto_dialog.setRange(area);

    auto all_area_xy = canvas_->getMapper().mapRange(canvas_->getDataRange());
    goto_dialog.setMinMaxOfRange(all_area_xy);
    // feature numbers only for consensus&feature maps
    goto_dialog.enableFeatureNumber(canvas()->getCurrentLayer().type == LayerDataBase::DT_FEATURE || canvas()->getCurrentLayer().type == LayerDataBase::DT_CONSENSUS);
    // execute
    if (goto_dialog.exec())
    {
      if (goto_dialog.showRange())
      {
        canvas()->setVisibleArea(goto_dialog.getRange());
      }
      else
      {
        String feature_id = goto_dialog.getFeatureNumber();
        //try to convert to UInt64 id
        UniqueIdInterface uid;
        uid.setUniqueId(feature_id);

        Size feature_index(-1); // TODO : not use -1
        auto* lf = dynamic_cast<LayerDataFeature*>(&canvas()->getCurrentLayer());
        auto* lc = dynamic_cast<LayerDataConsensus*>(&canvas()->getCurrentLayer());
        if (lf)
        {
          feature_index = lf->getFeatureMap()->uniqueIdToIndex(uid.getUniqueId());
        }
        else if (lc)
        {
          feature_index = lc->getConsensusMap()->uniqueIdToIndex(uid.getUniqueId());
        }
        if (feature_index == Size(-1)) // UID does not exist
        {
          try
          {
            feature_index = feature_id.toInt(); // normal feature index as stored in map
          }
          catch (...) // we might still deal with a UID, so toInt() will throw as the number is too big
          {
            feature_index = Size(-1);
          }
        }

        //check if the feature index exists
        if ((lf && feature_index >= lf->getFeatureMap()->size()) || (lc && feature_index >= lc->getConsensusMap()->size()))
        {
          QMessageBox::warning(this, "Invalid feature number", "Feature number too large/UniqueID not found.\nPlease select a valid feature!");
          return;
        }
        // display feature with a margin
        RangeAllType range;
        if (lf)
        {
          const FeatureMap& map = *lf->getFeatureMap();
          const DBoundingBox<2> bb = map[feature_index].getConvexHull().getBoundingBox();
          range.RangeRT::operator=(RangeBase{bb.minPosition()[0], bb.maxPosition()[0]});
          range.RangeMZ::operator=(RangeBase{bb.minPosition()[1], bb.maxPosition()[1]});
        }
        else // Consensus Feature
        {
          const ConsensusFeature& cf = (*lc->getConsensusMap())[feature_index];
          range = canvas_->getMapper().fromXY(canvas_->getMapper().map(cf));
        }
        range.RangeRT::extendLeftRight(30);
        range.RangeMZ::extendLeftRight(5);
        canvas()->setVisibleArea(range);
      }
    }
  }

  bool Plot2DWidget::projectionsVisible() const
  {
    return projection_onto_Y_->isVisible() || projection_onto_X_->isVisible();
  }

  void Plot2DWidget::autoUpdateProjections_()
  {
    if (projectionsVisible() && projections_auto_->isChecked())
    {
      projections_timer_->start();
    }
  }

} //Namespace
