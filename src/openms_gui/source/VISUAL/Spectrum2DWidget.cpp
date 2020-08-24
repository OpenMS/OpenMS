// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>

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

  Spectrum2DWidget::Spectrum2DWidget(const Param& preferences, QWidget* parent) :
    SpectrumWidget(preferences, parent)
  {
    setCanvas_(new Spectrum2DCanvas(preferences, this), 1, 2);

    // add axes
    x_axis_->setLegend(SpectrumWidget::MZ_AXIS_TITLE);
    y_axis_->setLegend(SpectrumWidget::RT_AXIS_TITLE);
    y_axis_->setMinimumWidth(50);

    // add projections
    grid_->setColumnStretch(2, 3);
    grid_->setRowStretch(1, 3);

    SpectrumCanvas::ExperimentSharedPtrType shr_ptr = SpectrumCanvas::ExperimentSharedPtrType(new SpectrumCanvas::ExperimentType());
    LayerData::ODExperimentSharedPtrType od_dummy(new OnDiscMSExperiment());
    MSSpectrum dummy_spec;
    dummy_spec.push_back(Peak1D());
    shr_ptr->addSpectrum(dummy_spec);

    projection_vert_ = new  Spectrum1DWidget(Param(), this);
    projection_vert_->hide();
    projection_vert_->canvas()->addLayer(shr_ptr, od_dummy);
    grid_->addWidget(projection_vert_, 1, 3, 2, 1);

    projection_horz_ = new Spectrum1DWidget(Param(), this);
    projection_horz_->canvas()->addLayer(shr_ptr, od_dummy);
    projection_horz_->hide();
    grid_->addWidget(projection_horz_, 0, 1, 1, 2);

    if (canvas()->isMzToXAxis())
    {
      projection_horz_->canvas()->setDrawMode(Spectrum1DCanvas::DM_PEAKS);
      projection_horz_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
      projection_vert_->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
      projection_vert_->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);
    }
    else
    {
      projection_horz_->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
      projection_horz_->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);
      projection_vert_->canvas()->setDrawMode(Spectrum1DCanvas::DM_PEAKS);
      projection_vert_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
    }

    connect(canvas(), SIGNAL(showProjectionHorizontal(ExperimentSharedPtrType)), this, SLOT(horizontalProjection(ExperimentSharedPtrType)));
    connect(canvas(), SIGNAL(showProjectionVertical(ExperimentSharedPtrType)), this, SLOT(verticalProjection(ExperimentSharedPtrType)));
    connect(canvas(), SIGNAL(showProjectionInfo(int, double, double)), this, SLOT(projectionInfo(int, double, double)));
    connect(canvas(), SIGNAL(toggleProjections()), this, SLOT(toggleProjections()));
    connect(canvas(), SIGNAL(visibleAreaChanged(DRange<2>)), this, SLOT(autoUpdateProjections()));
    // delegate signals from canvas
    connect(canvas(), SIGNAL(showSpectrumAs1D(int)), this, SIGNAL(showSpectrumAs1D(int)));
    connect(canvas(), SIGNAL(showSpectrumAs1D(std::vector<int, std::allocator<int> >)), this, SIGNAL(showSpectrumAs1D(std::vector<int, std::allocator<int> >)));
    connect(canvas(), SIGNAL(showCurrentPeaksAs3D()), this, SIGNAL(showCurrentPeaksAs3D()));
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
    connect(button, SIGNAL(clicked()), canvas(), SLOT(updateProjections()));
    box_grid->addWidget(button, 4, 0);

    projections_auto_ = new QCheckBox("Auto-update", projection_box_);
    projections_auto_->setWhatsThis("When activated, projections are automatically updated one second after the last change of the visible area.");
    projections_auto_->setChecked(true);
    box_grid->addWidget(projections_auto_, 4, 1);

    //set up projections auto-update
    projections_timer_ = new QTimer(this);
    projections_timer_->setSingleShot(true);
    projections_timer_->setInterval(1000);
    connect(projections_timer_, SIGNAL(timeout()), this, SLOT(updateProjections()));
  }

  Spectrum2DWidget::~Spectrum2DWidget()
  {
    //cout << "DEST Spectrum2DWidget" << endl;
  }

  void Spectrum2DWidget::projectionInfo(int peaks, double intensity, double max)
  {
    projection_peaks_->setText(QString::number(peaks));
    projection_sum_->setText(QString::number(intensity, 'f', 1));
    projection_max_->setText(QString::number(max, 'f', 1));
  }

  void Spectrum2DWidget::recalculateAxes_()
  {
    const SpectrumCanvas::AreaType area = canvas()->getVisibleArea();

    if (canvas()->isMzToXAxis())
    {
      x_axis_->setAxisBounds(area.minX(), area.maxX());
      y_axis_->setAxisBounds(area.minY(), area.maxY());
    }
    else
    {
      x_axis_->setAxisBounds(area.minY(), area.maxY());
      y_axis_->setAxisBounds(area.minX(), area.maxX());
    }
  }

  Histogram<> Spectrum2DWidget::createIntensityDistribution_() const
  {
    //initialize histogram
    double min = canvas_->getCurrentMinIntensity();
    double max = canvas_->getCurrentMaxIntensity();
    if (min == max)
    {
      min -= 0.01;
      max += 0.01;
    }
    Histogram<> tmp(min, max, (max - min) / 500.0);

    if (canvas_->getCurrentLayer().type == LayerData::DT_PEAK)
    {
      for (ExperimentType::ConstIterator spec_it = canvas_->getCurrentLayer().getPeakData()->begin(); spec_it != canvas_->getCurrentLayer().getPeakData()->end(); ++spec_it)
      {
        if (spec_it->getMSLevel() != 1)
        {
          continue;
        }
        for (ExperimentType::SpectrumType::ConstIterator peak_it = spec_it->begin(); peak_it != spec_it->end(); ++peak_it)
        {
          tmp.inc(peak_it->getIntensity());
        }
      }
    }
    else if (canvas_->getCurrentLayer().type == LayerData::DT_FEATURE)
    {
      for (Spectrum2DCanvas::FeatureMapType::ConstIterator it = canvas_->getCurrentLayer().getFeatureMap()->begin(); it != canvas_->getCurrentLayer().getFeatureMap()->end(); ++it)
      {
        tmp.inc(it->getIntensity());
      }
    }
    else
    {
      for (Spectrum2DCanvas::ConsensusMapType::ConstIterator it = canvas_->getCurrentLayer().getConsensusMap()->begin(); it != canvas_->getCurrentLayer().getConsensusMap()->end(); ++it)
      {
        tmp.inc(it->getIntensity());
      }
    }

    return tmp;
  }

  Histogram<> Spectrum2DWidget::createMetaDistribution_(const String& name) const
  {
    Histogram<> tmp;

    if (canvas_->getCurrentLayer().type == LayerData::DT_PEAK)
    {
      //determine min and max of the data
      float min = numeric_limits<float>::max(), max = -numeric_limits<float>::max();
      for (ExperimentType::const_iterator s_it = canvas_->getCurrentLayer().getPeakData()->begin(); s_it != canvas_->getCurrentLayer().getPeakData()->end(); ++s_it)
      {
        if (s_it->getMSLevel() != 1)
          continue;
        //float arrays
        for (ExperimentType::SpectrumType::FloatDataArrays::const_iterator it = s_it->getFloatDataArrays().begin(); it != s_it->getFloatDataArrays().end(); ++it)
        {
          if (it->getName() == name)
          {
            for (Size i = 0; i < it->size(); ++i)
            {
              if ((*it)[i] < min)
                min = (*it)[i];
              if ((*it)[i] > max)
                max = (*it)[i];
            }
            break;
          }
        }
        //integer arrays
        for (ExperimentType::SpectrumType::IntegerDataArrays::const_iterator it = s_it->getIntegerDataArrays().begin(); it != s_it->getIntegerDataArrays().end(); ++it)
        {
          if (it->getName() == name)
          {
            for (Size i = 0; i < it->size(); ++i)
            {
              if ((*it)[i] < min)
                min = (*it)[i];
              if ((*it)[i] > max)
                max = (*it)[i];
            }
            break;
          }
        }
      }
      if (min >= max)
        return tmp;

      //create histogram
      tmp.reset(min, max, (max - min) / 500.0);
      for (ExperimentType::const_iterator s_it = canvas_->getCurrentLayer().getPeakData()->begin(); s_it != canvas_->getCurrentLayer().getPeakData()->end(); ++s_it)
      {
        if (s_it->getMSLevel() != 1)
          continue;
        //float arrays
        for (ExperimentType::SpectrumType::FloatDataArrays::const_iterator it = s_it->getFloatDataArrays().begin(); it != s_it->getFloatDataArrays().end(); ++it)
        {
          if (it->getName() == name)
          {
            for (Size i = 0; i < it->size(); ++i)
            {
              tmp.inc((*it)[i]);
            }
            break;
          }
        }
        //integer arrays
        for (ExperimentType::SpectrumType::IntegerDataArrays::const_iterator it = s_it->getIntegerDataArrays().begin(); it != s_it->getIntegerDataArrays().end(); ++it)
        {
          if (it->getName() == name)
          {
            for (Size i = 0; i < it->size(); ++i)
            {
              tmp.inc((*it)[i]);
            }
            break;
          }
        }
      }
    }
    else //Features
    {
      //determine min and max
      float min = numeric_limits<float>::max(), max = -numeric_limits<float>::max();
      for (Spectrum2DCanvas::FeatureMapType::ConstIterator it = canvas_->getCurrentLayer().getFeatureMap()->begin(); it != canvas_->getCurrentLayer().getFeatureMap()->end(); ++it)
      {
        if (it->metaValueExists(name))
        {
          float value = it->getMetaValue(name);
          if (value < min)
            min = value;
          if (value > max)
            max = value;
        }
      }
      //create histogram
      tmp.reset(min, max, (max - min) / 500.0);
      for (Spectrum2DCanvas::FeatureMapType::ConstIterator it = canvas_->getCurrentLayer().getFeatureMap()->begin(); it != canvas_->getCurrentLayer().getFeatureMap()->end(); ++it)
      {
        if (it->metaValueExists(name))
        {
          tmp.inc((float)(it->getMetaValue(name)));
        }
      }

    }

    return tmp;
  }

  void Spectrum2DWidget::updateProjections()
  {
    canvas()->updateProjections();
  }

  void Spectrum2DWidget::toggleProjections()
  {
    if (projectionsVisible())
    {
      setMinimumSize(250, 250);
      projection_box_->hide();
      projection_horz_->hide();
      projection_vert_->hide();
      grid_->setColumnStretch(3, 0);
      grid_->setRowStretch(0, 0);
    }
    else
    {
      setMinimumSize(500, 500);
      updateProjections();
    }
  }

  void Spectrum2DWidget::horizontalProjection(ExperimentSharedPtrType exp)
  {
    LayerData::ODExperimentSharedPtrType od_dummy(new OnDiscMSExperiment());

    // print horizontal (note that m/z in the projection could actually be RT - this only determines the orientation)
    projection_horz_->canvas()->mzToXAxis(true); 
    projection_horz_->canvas()->setSwappedAxis(true);
    
    projection_horz_->showLegend(false);
    Spectrum1DCanvas::IntensityModes intensity = projection_horz_->canvas()->getIntensityMode();
    projection_horz_->canvas()->setIntensityMode(intensity);

    projection_horz_->canvas()->removeLayer(0);
    projection_horz_->canvas()->addLayer(exp, od_dummy);

    grid_->setColumnStretch(3, 2);

        if (canvas()->isMzToXAxis())
    {
      projection_horz_->canvas()->setDrawMode(Spectrum1DCanvas::DM_PEAKS);
      projection_horz_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
      projection_vert_->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
      projection_vert_->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);
    }
    else
    {
      projection_horz_->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
      projection_horz_->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);
      projection_vert_->canvas()->setDrawMode(Spectrum1DCanvas::DM_PEAKS);
      projection_vert_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
    }
    projection_horz_->show();
    projection_box_->show();
  }

  void Spectrum2DWidget::verticalProjection(ExperimentSharedPtrType exp)
  {
    LayerData::ODExperimentSharedPtrType od_dummy(new OnDiscMSExperiment());
    // print vertically (note that m/z in the projection could actually be RT - this only determines the orientation)
    projection_vert_->canvas()->mzToXAxis(false); 
    projection_vert_->canvas()->setSwappedAxis(true);

    projection_vert_->showLegend(false);
    Spectrum1DCanvas::IntensityModes intensity = projection_vert_->canvas()->getIntensityMode();
    projection_vert_->canvas()->setIntensityMode(intensity);

    projection_vert_->canvas()->removeLayer(0);
    projection_vert_->canvas()->addLayer(exp, od_dummy);

    grid_->setRowStretch(0, 2);

    if (canvas()->isMzToXAxis())
    {
      projection_horz_->canvas()->setDrawMode(Spectrum1DCanvas::DM_PEAKS);
      projection_horz_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
      projection_vert_->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
      projection_vert_->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);
    }
    else
    {
      projection_horz_->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
      projection_horz_->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);
      projection_vert_->canvas()->setDrawMode(Spectrum1DCanvas::DM_PEAKS);
      projection_vert_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
    }
    projection_box_->show();
    projection_vert_->show();
  }

  const Spectrum1DWidget* Spectrum2DWidget::getHorizontalProjection() const
  {
    return projection_horz_;
  }

  const Spectrum1DWidget* Spectrum2DWidget::getVerticalProjection() const
  {
    return projection_vert_;
  }

  void Spectrum2DWidget::showGoToDialog()
  {
    Spectrum2DGoToDialog goto_dialog(this);
    //set range
    const DRange<2>& area = canvas()->getVisibleArea();
    goto_dialog.setRange(area.minY(), area.maxY(), area.minX(), area.maxX());
    goto_dialog.setMinMaxOfRange(canvas()->getDataRange().minY(), canvas()->getDataRange().maxY(), canvas()->getDataRange().minX(), canvas()->getDataRange().maxX());
    // feature numbers only for consensus&feature maps
    goto_dialog.enableFeatureNumber(canvas()->getCurrentLayer().type == LayerData::DT_FEATURE || canvas()->getCurrentLayer().type == LayerData::DT_CONSENSUS);
    //execute
    if (goto_dialog.exec())
    {
      if (goto_dialog.showRange())
      {
        goto_dialog.fixRange();
        SpectrumCanvas::AreaType area(goto_dialog.getMinMZ(), goto_dialog.getMinRT(), goto_dialog.getMaxMZ(), goto_dialog.getMaxRT());
        if (goto_dialog.checked()) correctAreaToObeyMinMaxRanges_(area);
        canvas()->setVisibleArea(area);
      }
      else
      {
        String feature_id = goto_dialog.getFeatureNumber();
        //try to convert to UInt64 id
        UniqueIdInterface uid;
        uid.setUniqueId(feature_id);

        Size feature_index(-1); // TODO : not use -1
        if (canvas()->getCurrentLayer().type == LayerData::DT_FEATURE)
          feature_index = canvas()->getCurrentLayer().getFeatureMap()->uniqueIdToIndex(uid.getUniqueId());
        else if (canvas()->getCurrentLayer().type == LayerData::DT_CONSENSUS)
          feature_index = canvas()->getCurrentLayer().getConsensusMap()->uniqueIdToIndex(uid.getUniqueId());
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
        if ((canvas()->getCurrentLayer().type == LayerData::DT_FEATURE && feature_index >= canvas()->getCurrentLayer().getFeatureMap()->size())
           || (canvas()->getCurrentLayer().type == LayerData::DT_CONSENSUS && feature_index >= canvas()->getCurrentLayer().getConsensusMap()->size()))
        {
          QMessageBox::warning(this, "Invalid feature number", "Feature number too large/UniqueID not found.\nPlease select a valid feature!");
          return;
        }
        //display feature with a margin
        if (canvas()->getCurrentLayer().type == LayerData::DT_FEATURE)
        {
          const FeatureMapType& map = *canvas()->getCurrentLayer().getFeatureMap();
          DBoundingBox<2> bb = map[feature_index].getConvexHull().getBoundingBox();
          double rt_margin = (bb.maxPosition()[0] - bb.minPosition()[0]) * 0.5;
          double mz_margin = (bb.maxPosition()[1] - bb.minPosition()[1]) * 2;
          SpectrumCanvas::AreaType narea(bb.minPosition()[1] - mz_margin, bb.minPosition()[0] - rt_margin, bb.maxPosition()[1] + mz_margin, bb.maxPosition()[0] + rt_margin);
          canvas()->setVisibleArea(narea);
        }
        else // Consensus Feature
        {
          const ConsensusFeature& cf = (*canvas()->getCurrentLayer().getConsensusMap())[feature_index];
          double rt_margin = 30;
          double mz_margin = 5;
          SpectrumCanvas::AreaType narea(cf.getMZ() - mz_margin, cf.getRT() - rt_margin, cf.getMZ() + mz_margin, cf.getRT() + rt_margin);
          canvas()->setVisibleArea(narea);
        }

      }
    }
  }

  bool Spectrum2DWidget::projectionsVisible() const
  {
    return projection_horz_->isVisible() || projection_vert_->isVisible();
  }

  void Spectrum2DWidget::autoUpdateProjections()
  {
    if (projectionsVisible() && projections_auto_->isChecked())
    {
      projections_timer_->start();
    }
  }

} //Namespace
