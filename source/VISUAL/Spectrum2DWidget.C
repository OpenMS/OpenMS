// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>

#include <QtGui/QPushButton>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QMessageBox>
#include <QtGui/QCheckBox>
#include <QtCore/QTimer>

using namespace std;

namespace OpenMS
{
	using namespace Internal;
	using namespace Math;
	
	Spectrum2DWidget::Spectrum2DWidget(const Param& preferences, QWidget* parent)
		: SpectrumWidget(preferences, parent)
	{
		setCanvas_(new Spectrum2DCanvas(preferences, this),1,2);
		
		// add axes
		x_axis_->setLegend(String(Peak2D::shortDimensionName(Peak2D::MZ))+" ["+String(Peak2D::shortDimensionUnit(Peak2D::MZ))+"]");
		y_axis_->setLegend(String(Peak2D::shortDimensionName(Peak2D::RT))+" ["+String(Peak2D::shortDimensionUnit(Peak2D::RT))+"]");
		y_axis_->setMinimumWidth(50);
		
		// add projetions
		grid_->setColumnStretch(2,3);
		grid_->setRowStretch(1,3);
		
		projection_vert_ = new 	Spectrum1DWidget(Param(), this);
		projection_vert_->hide();
		grid_->addWidget(projection_vert_,1,3,2,1);
		
		projection_horz_ = new Spectrum1DWidget(Param(), this);
		projection_horz_->hide();
		grid_->addWidget(projection_horz_,0,1,1,2);
    connect(canvas(), SIGNAL(showProjectionHorizontal(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes)), this, SLOT(horizontalProjection(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes)));
    connect(canvas(), SIGNAL(showProjectionVertical(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes)), this, SLOT(verticalProjection(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes)));
		connect(canvas(), SIGNAL(showProjectionInfo(int,double,double)), this, SLOT(projectionInfo(int,double,double)));
		connect(canvas(), SIGNAL(toggleProjections()), this, SLOT(toggleProjections()));
		connect(canvas(), SIGNAL(visibleAreaChanged(DRange<2>)), this, SLOT(autoUpdateProjections()));
    // delegate signals from canvas
    connect(canvas(), SIGNAL(showSpectrumAs1D(int)), this, SIGNAL(showSpectrumAs1D(int)));
    connect(canvas(), SIGNAL(showCurrentPeaksAs3D()), this, SIGNAL(showCurrentPeaksAs3D()));
		// add projections box
		projection_box_ = new QGroupBox("Projections",this);
		projection_box_->hide();
		grid_->addWidget(projection_box_, 0, 3);
		QGridLayout* box_grid = new QGridLayout(projection_box_); 
		
		QLabel* label = new QLabel("Peaks: ");
		box_grid->addWidget(label,0,0);
		projection_peaks_ = new QLabel("");
		box_grid->addWidget(projection_peaks_,0,1);
		
		label = new QLabel("Intensity sum: ");
		box_grid->addWidget(label,1,0);
		projection_sum_ = new QLabel("");
		box_grid->addWidget(projection_sum_,1,1);

		label = new QLabel("Maximum intensity: ");
		box_grid->addWidget(label,2,0);
		projection_max_ = new QLabel("");
		box_grid->addWidget(projection_max_,2,1);
		
		box_grid->setRowStretch(3,2);

		QPushButton* button = new QPushButton("Update", projection_box_);
		connect(button, SIGNAL(clicked()), canvas(), SLOT(updateProjections()));
		box_grid->addWidget(button,4,0);

		projections_auto_ = new QCheckBox("Auto-update", projection_box_);
		projections_auto_->setWhatsThis("When activated, projections are automatically updated one second after the last change of the visible area.");
		projections_auto_->setChecked(true);
		box_grid->addWidget(projections_auto_,4,1);
		
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
		projection_sum_->setText(QString::number(intensity,'f',1));
		projection_max_->setText(QString::number(max,'f',1));
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
		DoubleReal min = canvas_->getCurrentMinIntensity();
		DoubleReal max = canvas_->getCurrentMaxIntensity();
		if (min==max)
		{
			min-=0.01;
			max+=0.01;
		}
		Histogram<> tmp(min,max,(max-min)/500.0);
		
		if (canvas_->getCurrentLayer().type==LayerData::DT_PEAK)
		{
      for (ExperimentType::ConstIterator spec_it = canvas_->getCurrentLayer().getPeakData()->begin(); spec_it != canvas_->getCurrentLayer().getPeakData()->end(); ++spec_it)
			{
				if (spec_it->getMSLevel()!=1)
				{
					continue;
				}
				for (ExperimentType::SpectrumType::ConstIterator peak_it = spec_it->begin(); peak_it != spec_it->end(); ++peak_it)
				{
					tmp.inc(peak_it->getIntensity());
				}
			}
		}
		else if (canvas_->getCurrentLayer().type==LayerData::DT_FEATURE)
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

		if (canvas_->getCurrentLayer().type==LayerData::DT_PEAK)
		{
			//determine min and max of the data
			Real min = numeric_limits<Real>::max(), max = -numeric_limits<Real>::max();
                        for (ExperimentType::const_iterator s_it = canvas_->getCurrentLayer().getPeakData()->begin(); s_it!=canvas_->getCurrentLayer().getPeakData()->end(); ++s_it)
			{
				if (s_it->getMSLevel()!=1) continue;
				//float arrays
				for (ExperimentType::SpectrumType::FloatDataArrays::const_iterator it=s_it->getFloatDataArrays().begin(); it!=s_it->getFloatDataArrays().end(); it++)
				{
					if (it->getName()==name)
					{
						for (Size i=0; i<it->size(); ++i)
						{
							if ((*it)[i]<min) min = (*it)[i];
							if ((*it)[i]>max) max = (*it)[i];
						}
						break;
					}
				}
				//integer arrays
				for (ExperimentType::SpectrumType::IntegerDataArrays::const_iterator it=s_it->getIntegerDataArrays().begin(); it!=s_it->getIntegerDataArrays().end(); it++)
				{
					if (it->getName()==name)
					{
						for (Size i=0; i<it->size(); ++i)
						{
							if ((*it)[i]<min) min = (*it)[i];
							if ((*it)[i]>max) max = (*it)[i];
						}
						break;
					}
				}
			}
			if (min>=max) return tmp;
		
			//create histogram
			tmp.reset(min,max,(max-min)/500.0);
                        for (ExperimentType::const_iterator s_it = canvas_->getCurrentLayer().getPeakData()->begin(); s_it!=canvas_->getCurrentLayer().getPeakData()->end(); ++s_it)
			{
				if (s_it->getMSLevel()!=1) continue;
				//float arrays
				for (ExperimentType::SpectrumType::FloatDataArrays::const_iterator it=s_it->getFloatDataArrays().begin(); it!=s_it->getFloatDataArrays().end(); it++)
				{
					if (it->getName()==name)
					{
						for (Size i=0; i<it->size(); ++i)
						{
							tmp.inc((*it)[i]);
						}
						break;
					}
				}
				//integer arrays
				for (ExperimentType::SpectrumType::IntegerDataArrays::const_iterator it=s_it->getIntegerDataArrays().begin(); it!=s_it->getIntegerDataArrays().end(); it++)
				{
					if (it->getName()==name)
					{
						for (Size i=0; i<it->size(); ++i)
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
			Real min = numeric_limits<Real>::max(), max = -numeric_limits<Real>::max();
      for (Spectrum2DCanvas::FeatureMapType::ConstIterator it = canvas_->getCurrentLayer().getFeatureMap()->begin(); it != canvas_->getCurrentLayer().getFeatureMap()->end(); ++it)
			{
				if (it->metaValueExists(name))
				{
					Real value = it->getMetaValue(name);
					if (value<min) min = value;
					if (value>max) max = value;
				}
			}
			//create histogram
			tmp.reset(min,max,(max-min)/500.0);
      for (Spectrum2DCanvas::FeatureMapType::ConstIterator it = canvas_->getCurrentLayer().getFeatureMap()->begin(); it != canvas_->getCurrentLayer().getFeatureMap()->end(); ++it)
			{
				if (it->metaValueExists(name))
				{
					tmp.inc((Real)(it->getMetaValue(name)));
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
			setMinimumSize(250,250);
			projection_box_->hide();
			projection_horz_->hide();
			projection_vert_->hide();
			grid_->setColumnStretch(3,0);
			grid_->setRowStretch(0,0);
		}
		else
		{
			setMinimumSize(500,500);
			updateProjections();
		}
	}
	
  void Spectrum2DWidget::horizontalProjection(ExperimentSharedPtrType exp, Spectrum1DCanvas::DrawModes mode)
	{
		projection_horz_->canvas()->mzToXAxis(true); // determines the orientation of the data
    projection_horz_->canvas()->setSwappedAxis(canvas()->isMzToXAxis());
		projection_horz_->showLegend(false);
		projection_horz_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
		if (!projectionsVisible() || projection_horz_->canvas()->getLayerCount()==0) //set draw mode
		{
			projection_horz_->canvas()->removeLayer(0);
      projection_horz_->canvas()->addLayer(exp);
			projection_horz_->canvas()->setDrawMode(mode);
		}
		else //keep draw mode
		{
			Spectrum1DCanvas::DrawModes previous_mode = projection_horz_->canvas()->getDrawMode();
			projection_horz_->canvas()->removeLayer(0);
      projection_horz_->canvas()->addLayer(exp);
			projection_horz_->canvas()->setDrawMode(previous_mode);
		}
		grid_->setColumnStretch(3,2);
		projection_horz_->show();
		projection_box_->show();
	}
	
  void Spectrum2DWidget::verticalProjection(ExperimentSharedPtrType exp, Spectrum1DCanvas::DrawModes mode)
	{
		projection_vert_->canvas()->mzToXAxis(false); // determines the orientation of the data
    projection_vert_->canvas()->setSwappedAxis(canvas()->isMzToXAxis());
		projection_vert_->showLegend(false);
		projection_vert_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
    // todo: why would we want to retain the old draw mode (same in horizontalProjection)??
		if (!projectionsVisible() || projection_vert_->canvas()->getLayerCount()==0) //set draw mode
		{
			projection_vert_->canvas()->removeLayer(0);
      projection_vert_->canvas()->addLayer(exp);
			projection_vert_->canvas()->setDrawMode(mode);
		}
		else //keep draw mode
		{
			Spectrum1DCanvas::DrawModes previous_mode = projection_vert_->canvas()->getDrawMode();
			projection_vert_->canvas()->removeLayer(0);
      projection_vert_->canvas()->addLayer(exp);
			projection_vert_->canvas()->setDrawMode(previous_mode);
		}
		grid_->setRowStretch(0,2);
		projection_vert_->show();
		projection_box_->show();
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
	  goto_dialog.setRange(area.minY(),area.maxY(),area.minX(),area.maxX()); 
	  // feature numbers only for consensus&feature maps
	  goto_dialog.enableFeatureNumber(canvas()->getCurrentLayer().type == LayerData::DT_FEATURE || canvas()->getCurrentLayer().type == LayerData::DT_CONSENSUS);
	  //execute
	  if (goto_dialog.exec())
	  {
	  	if (goto_dialog.showRange())
	  	{
	  		canvas()->setVisibleArea(SpectrumCanvas::AreaType( goto_dialog.getMinMZ(), goto_dialog.getMinRT(), goto_dialog.getMaxMZ(), goto_dialog.getMaxRT()));
			}
			else
			{
        String feature_id = goto_dialog.getFeatureNumber();
        //try to convert to UInt64 id
        UniqueIdInterface uid;
        uid.setUniqueId(feature_id);

        Size feature_index;
        if (canvas()->getCurrentLayer().type == LayerData::DT_FEATURE) feature_index = canvas()->getCurrentLayer().getFeatureMap()->uniqueIdToIndex(uid.getUniqueId());
        else if (canvas()->getCurrentLayer().type == LayerData::DT_CONSENSUS)  feature_index = canvas()->getCurrentLayer().getConsensusMap()->uniqueIdToIndex(uid.getUniqueId());
        if (feature_index == Size(-1)) // UID does not exist
        {
          try
          {
            feature_index=feature_id.toInt(); // normal feature index as stored in map
          }
          catch (...)
          { // we might still deal with a UID, so toInt() will throw as the number is too big
            feature_index = Size(-1);
          }
        }

				//check if the feature index exists
        if ( (canvas()->getCurrentLayer().type == LayerData::DT_FEATURE && feature_index>=canvas()->getCurrentLayer().getFeatureMap()->size())
           ||(canvas()->getCurrentLayer().type == LayerData::DT_CONSENSUS && feature_index>=canvas()->getCurrentLayer().getConsensusMap()->size()))
				{
					QMessageBox::warning(this, "Invalid feature number", "Feature number too large/UniqueID not found.\nPlease select a valid feature!");
					return;
				}
				//display feature with a margin
        if (canvas()->getCurrentLayer().type == LayerData::DT_FEATURE)
        {
          const FeatureMapType& map = *canvas()->getCurrentLayer().getFeatureMap();
          DBoundingBox<2> bb = map[feature_index].getConvexHull().getBoundingBox();
          DoubleReal rt_margin = (bb.maxPosition()[0] - bb.minPosition()[0])*0.5;
          DoubleReal mz_margin = (bb.maxPosition()[1] - bb.minPosition()[1])*2;
          SpectrumCanvas::AreaType narea(bb.minPosition()[1]-mz_margin, bb.minPosition()[0]-rt_margin, bb.maxPosition()[1]+mz_margin, bb.maxPosition()[0]+rt_margin);
          canvas()->setVisibleArea(narea);
        }
        else // Consensus Feature
        {
          const ConsensusFeature& cf = (*canvas()->getCurrentLayer().getConsensusMap())[feature_index];
          DoubleReal rt_margin = 30;
          DoubleReal mz_margin = 5;
          SpectrumCanvas::AreaType narea(cf.getMZ()-mz_margin, cf.getRT()-rt_margin, cf.getMZ()+mz_margin, cf.getRT()+rt_margin);
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

