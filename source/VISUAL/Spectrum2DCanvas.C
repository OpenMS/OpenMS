// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/CONCEPT/TimeStamp.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DPrefDialog.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>

//STL
#include <algorithm>	

//QT
#include <QtGui/QWheelEvent>
#include <QtGui/QMouseEvent>
#include <QtGui/QPainter>
#include <QtGui/QMenu>
#include <QtGui/QBitmap>
#include <QtGui/QPolygon>
#include <QtCore/QTime>
#include <QtGui/QComboBox>


using namespace std;

namespace OpenMS
{
	using namespace Internal;

	Spectrum2DCanvas::Spectrum2DCanvas(const Param& preferences, QWidget* parent)
		: SpectrumCanvas(preferences, parent),
			selected_peak_(0),
			measurement_start_(0),
			measurement_stop_(0),
			tmp_peak_()
	{
    //Paramater handling
    defaults_.setValue("background_color", "#ffffff", "Background color.");
    defaults_.setValue("interpolation_steps", 200, "Number of interploation steps for peak gradient precalculation.");
    defaults_.setValue("dot:gradient", "Linear|0,#efef00;7,#ffaa00;15,#ff0000;27,#aa00ff;55,#5500ff;100,#000000", "Multi-color gradient for peaks.");
    defaults_.setValue("mapping_of_mz_to","x_axis","Determines with axis is the m/z axis.");
		defaultsToParam_();
		setName("Spectrum2DCanvas");
		setParameters(preferences);

		projection_mz_.resize(1);
		projection_rt_.resize(1);
		
		//set preferences and update widgets acoordningly
		if (param_.getValue("mapping_of_mz_to") != "x_axis")
		{
			mzToXAxis(false);
		}
	}
	
	Spectrum2DCanvas::~Spectrum2DCanvas()
	{
		//cout << "DEST Spectrum2DCanvas" << endl;
	}
	
	void Spectrum2DCanvas::highlightPeak_(QPainter& painter, Feature* peak)
	{
		if (!peak) return;
		painter.save();
		painter.setPen(QPen(Qt::red, 2));
		QPoint pos;
		dataToWidget_(peak->getMZ(),peak->getRT(),pos);
		painter.drawEllipse(pos.x() - 5, pos.y() - 5, 10, 10);
		painter.restore();
	}
	
	Feature* Spectrum2DCanvas::findNearestPeak_(const QPoint& pos)
	{
		//Constructing the area corrects swapped mapping of RT and m/z
		AreaType area (widgetToData_(pos - QPoint(5,5)),widgetToData_(pos + QPoint(5,5)));

		Feature* max_peak = 0;
		float max_int = -1 * numeric_limits<float>::max();
		
		//cout << "findNearestPeak_: Int range -- " << getCurrentLayer().min_int << " "  << getCurrentLayer().max_int << endl;
		
		if (getCurrentLayer().type==LayerData::DT_PEAK)
		{
			for (ExperimentType::ConstAreaIterator i = getCurrentPeakData().areaBeginConst(area.min()[1],area.max()[1],area.min()[0],area.max()[0]); 
					 i != getCurrentPeakData().areaEndConst(); 
					 ++i)
			{
				if (i->getIntensity() > max_int && i->getIntensity()>=getCurrentLayer().min_int && i->getIntensity()<=getCurrentLayer().max_int)
				{
					//cout << "new max: " << i.getRT() << " " << i->getMZ() << endl;
					max_int = i->getIntensity();
					
					tmp_peak_.setIntensity(i->getIntensity());
					tmp_peak_.setMZ(i->getMZ());
					tmp_peak_.setRT(i.getRT());
					tmp_peak_.setCharge(0);
					max_peak = &tmp_peak_;
				}
			}
	 	}
	 	else
	 	{
			for (FeatureMapType::ConstIterator i = getCurrentLayer().features.begin();
				   i != getCurrentLayer().features.end();
				   ++i)
			{
				if ( i->getRT() >= area.min()[1] &&
						 i->getRT() <= area.max()[1] &&
						 i->getMZ() >= area.min()[0] &&
						 i->getMZ() <= area.max()[0] &&
						 i->getIntensity()    >= getCurrentLayer().min_int &&
						 i->getIntensity()    <= getCurrentLayer().max_int )
				{
					if (i->getIntensity() > max_int)
					{
						max_int = i->getIntensity();
						
						tmp_peak_.setIntensity(i->getIntensity());
						tmp_peak_.setMZ(i->getMZ());
						tmp_peak_.setRT(i->getRT());
						tmp_peak_.setCharge(i->getCharge());
						tmp_peak_.getConvexHulls() = i->getConvexHulls();
						
						max_peak = &tmp_peak_;
					}				
				}
			}	 	
		}
//		if (max_peak!=0)
//		{
//			cout << "MAX PEAK: " << max_peak->getMZ() << endl;
//		}
		return max_peak;
	}
	
	float Spectrum2DCanvas::betweenFactor_(float v1, float v2, float val)
	{
		float lo = min(v1, v2);
		float hi = max(v1, v2);
		return (hi - lo == 0) ? 1 : (val - lo) / (hi - lo);
	}
	
	const QColor& Spectrum2DCanvas::heightColor_(float val, const MultiGradient& gradient)
	{
		switch (intensity_mode_)
		{
			case IM_NONE:
				return gradient.precalculatedColorAt(val);
				break;
			case IM_LOG:
				return gradient.precalculatedColorAt(log(val+1)); //prevent log of numbers samller than 1
				break;
			case IM_PERCENTAGE:
				return gradient.precalculatedColorAt(val*percentage_factor_);
				break;
			case IM_SNAP:
				return gradient.precalculatedColorAt(val*snap_factor_);
				break;
			default:
				throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
	}

	void Spectrum2DCanvas::paintDots_(UInt layer_index, QPainter& painter)
	{
#ifdef TIMING_TOPPVIEW
		QTime timer;
		timer.start();
#endif	
		
		if (intensity_mode_ == IM_PERCENTAGE)
		{
			if (getLayer(layer_index).type == LayerData::DT_PEAK)
			{
				percentage_factor_ = overall_data_range_.max()[2]/getPeakData(layer_index).getMaxInt();
			}
			else
			{
				percentage_factor_ = overall_data_range_.max()[2]/getLayer(layer_index).features.getMaxInt();
			}
		}
		else 
		{
			percentage_factor_ = 1.0;
		}
		
		//temporary variable
		QPoint pos;
		
		double min_int = getLayer(layer_index).min_int;
		double max_int = getLayer(layer_index).max_int;
		
		painter.setPen(Qt::black);
		
		if (getLayer(layer_index).type==LayerData::DT_PEAK) //peaks
		{
			//determine number of MS1 scans
			UInt scans = 0;
			PeakMap::ConstIterator it = getPeakData(layer_index).RTBegin(visible_area_.min()[1]);
			while (it != getPeakData(layer_index).RTEnd(visible_area_.max()[1]))
			{
				if (it->getMSLevel()==1) ++scans;
				++it;
			}
			//determine number of shown peaks
			Int peaks = 0;
			it = getPeakData(layer_index).RTBegin(visible_area_.min()[1]) + scans/2;
			while (it!=getPeakData(layer_index).end() && it->getMSLevel()!=1)
			{
				++it;
			}
			if (it!=getPeakData(layer_index).end())
			{
				PeakSpectrum::ConstIterator it2 = it->MZBegin(visible_area_.min()[0]);
				while (it2!=it->MZEnd(visible_area_.max()[0]))
				{
					if (it2->getIntensity()>=min_int && it2->getIntensity()<=max_int) ++peaks;
					++it2;
				}
			}
			//paint dots if too many peaks are shown, crosses otherwise
			bool dots = false;
			if(isMzToXAxis())
			{
				if (peaks>0.5*width() || scans>0.5*height()) dots=true;
			}
			else
			{
				if (peaks>0.5*height() || scans>0.5*width()) dots=true;
			}
			//cout << "peaks: " << peaks << "  scans: " << scans << endl;
			//cout << "width: " << width() << "  height: " << height() << endl;
			for (ExperimentType::ConstAreaIterator i = getPeakData(layer_index).areaBeginConst(visible_area_.min()[1],visible_area_.max()[1],visible_area_.min()[0],visible_area_.max()[0]); 
					 i != getPeakData(layer_index).areaEndConst(); 
					 ++i)
			{
				if (i->getIntensity()>=min_int && i->getIntensity()<=max_int)
				{
					painter.setPen(heightColor_(i->getIntensity(), getLayer(layer_index).gradient));
					if (dots)
					{
						dataToWidget_(i->getMZ(), i.getRT(),pos);
						painter.drawPoint(pos.x(),pos.y());
					}
					else
					{
						dataToWidget_(i->getMZ(), i.getRT(),pos);
						painter.drawLine(pos.x(),pos.y()-1,pos.x(),pos.y()+1);
						painter.drawLine(pos.x()-1,pos.y(),pos.x()+1,pos.y());
					}
				}
			}
			
			//draw precursor peaks
			if (getLayerFlag(layer_index,LayerData::P_PRECURSORS))
			{
				const ExperimentType& exp = getPeakData(layer_index); 
				painter.setPen(Qt::black);
				for (ExperimentType::ConstIterator i = exp.RTBegin(visible_area_.min()[1]); 
						 i != exp.RTEnd(visible_area_.max()[1]); 
						 ++i)
				{
					//this is an MS/MS scan
					if (i->getMSLevel()==2)
					{
						ExperimentType::ConstIterator prec=exp.getPrecursorSpectrum(i);
						if (prec!=exp.end())
						{
							dataToWidget_(i->getPrecursorPeak().getPosition()[0], prec->getRT(),pos);
							painter.drawLine(pos.x(),pos.y()+2,pos.x()+2,pos.y());
							painter.drawLine(pos.x()+2,pos.y(),pos.x(),pos.y()-2);
							painter.drawLine(pos.x(),pos.y()-2,pos.x()-2,pos.y());
							painter.drawLine(pos.x()-2,pos.y(),pos.x(),pos.y()+2);
						}
					}
				}
			}
		}
		else //features
		{
			bool numbers = getLayerFlag(layer_index,LayerData::F_NUMBERS);
			UInt num=0;
			for (FeatureMapType::ConstIterator i = getLayer(layer_index).features.begin();
				   i != getLayer(layer_index).features.end();
				   ++i)
			{
				++num;
				if ( i->getRT() >= visible_area_.min()[1] &&
						 i->getRT() <= visible_area_.max()[1] &&
						 i->getMZ() >= visible_area_.min()[0] &&
						 i->getMZ() <= visible_area_.max()[0] &&
						 i->getIntensity()>=min_int &&
						 i->getIntensity()<=max_int)
				{
					painter.setPen(heightColor_(i->getIntensity(), getLayer(layer_index).gradient));
					dataToWidget_(i->getMZ(),i->getRT(),pos);
					painter.drawLine(pos.x(),pos.y()-1,pos.x(),pos.y()+1);
					painter.drawLine(pos.x()-1,pos.y(),pos.x()+1,pos.y());
					if (numbers)
					{
						painter.setPen(Qt::black);
						painter.drawText(pos.x()+10,pos.y()+10,QString::number(num));
					}
				}
			}
		}

#ifdef TIMING_TOPPVIEW
		cout << "paintDots_ took " << timer.elapsed() << " ms" << endl;
#endif	
	}

	void Spectrum2DCanvas::paintConvexHulls_(UInt layer_index, QPainter& painter)
	{
		painter.setPen(Qt::black);

		double min_int = getLayer(layer_index).min_int;
		double max_int = getLayer(layer_index).max_int;

		for (FeatureMapType::ConstIterator i = getLayer(layer_index).features.begin();
			   i != getLayer(layer_index).features.end();
			   ++i)
		{
			if ( i->getRT() >= visible_area_.min()[1] &&
					 i->getRT() <= visible_area_.max()[1] &&
					 i->getMZ() >= visible_area_.min()[0] &&
					 i->getMZ() <= visible_area_.max()[0] &&
					 i->getIntensity()>=min_int &&
					 i->getIntensity()<=max_int)
			{
				paintConvexHulls_(i->getConvexHulls(),painter);
			}
		}
	}            

	void Spectrum2DCanvas::paintFeaturePairConnections_(UInt layer_index, QPainter& painter)
	{
		painter.setPen(Qt::black);

		double min_int = getLayer(layer_index).min_int;
		double max_int = getLayer(layer_index).max_int;

		QPoint line_begin, line_end;
		FeatureMapType::ConstIterator i2;
	
		for (FeatureMapType::ConstIterator i1 = getLayer(layer_index).features.begin();
			   i1 != getLayer(layer_index).features.end();
			   i1+=2)
		{
			//get second feature
			i2 = i1 + 1;
			
			if ( i1->getRT() >= visible_area_.min()[1] &&
					 i1->getRT() <= visible_area_.max()[1] &&
					 i1->getMZ() >= visible_area_.min()[0] &&
					 i1->getMZ() <= visible_area_.max()[0] &&
					 i1->getIntensity()>=min_int &&
					 i1->getIntensity()<=max_int &&
					 i2->getRT() >= visible_area_.min()[1] &&
					 i2->getRT() <= visible_area_.max()[1] &&
					 i2->getMZ() >= visible_area_.min()[0] &&
					 i2->getMZ() <= visible_area_.max()[0] &&
					 i2->getIntensity()>=min_int &&
					 i2->getIntensity()<=max_int
					 )
			{
				dataToWidget_(i1->getMZ(),i1->getRT(), line_begin);
				dataToWidget_(i2->getMZ(),i2->getRT(), line_end);
				painter.drawLine(line_begin, line_end);
			}
		}
	}      

  void Spectrum2DCanvas::paintConvexHulls_(const Feature::ConvexHullVector& hulls, QPainter& painter)
  {
		QPolygon points;
		
		//iterate over all convex hulls
		for (UInt hull=0; hull<hulls.size(); ++hull)
		{
			points.resize(hulls[hull].getPoints().size());
			UInt index=0;
			QPoint pos;
			//iterate over hull points
			for(ConvexHull2D::PointArrayType::const_iterator it=hulls[hull].getPoints().begin(); it!=hulls[hull].getPoints().end(); ++it, ++index)
			{
				dataToWidget_(it->getY(), it->getX(),pos);
				points.setPoint(index, pos);
			}	
			//cout << "Hull: " << hull << " Points: " << points.size()<<endl;
			painter.drawPolygon(points);
		}
  }
	
	void Spectrum2DCanvas::intensityDistributionChange_()
	{
		update_buffer_ = true;
		update();
	}
	
	void Spectrum2DCanvas::intensityModeChange_()
	{
		for (UInt i=0; i<layers_.size();++i)
		{
			recalculateDotGradient_(i);
		}
		SpectrumCanvas::intensityModeChange_();
	}
	
	void Spectrum2DCanvas::recalculateDotGradient_(UInt layer)
	{
		getLayer_(layer).gradient.fromString(getLayer_(layer).param.getValue("dot:gradient"));
		//cout << "recalculateDotGradient_" << endl;
		if (intensity_mode_ == IM_LOG)
		{
			//cout << "LOG:" <<" "<< log(overall_data_range_.min()[2]) <<" "<< log(overall_data_range_.max()[2])<<" "<<param_.getValue("interpolation_steps")<<endl;
			getLayer_(layer).gradient.activatePrecalculationMode(0, log(overall_data_range_.max()[2]+1), param_.getValue("interpolation_steps"));
		}
		else
		{
			//cout << "NORMAL:" << overall_data_range_.min()[2] <<" "<< overall_data_range_.max()[2]<<" "<<param_.getValue("interpolation_steps")<<endl;
			getLayer_(layer).gradient.activatePrecalculationMode(0, overall_data_range_.max()[2], param_.getValue("interpolation_steps"));
		}	
	}
	
	void Spectrum2DCanvas::showProjections()
	{
		//create projection data
		map<float, float> rt;
		map<int, float> mzint;
		map<int,int> mzcount;
		map<int,float> mzsum;
    
		UInt peak_count = 0;
		DoubleReal intensity_max = 0.0;
		DoubleReal intensity_sum = 0.0;
	
	  float prec = 0.05;
		float mult = 1.0/prec;


		for (ExperimentType::ConstAreaIterator i = getCurrentPeakData().areaBeginConst(visible_area_.min()[1],visible_area_.max()[1],visible_area_.min()[0],visible_area_.max()[0]); 
				 i != getCurrentPeakData().areaEndConst();
				 ++i)
		{
			if (i->getIntensity()>=getCurrentLayer().min_int && i->getIntensity()<=getCurrentLayer().max_int)
			{
				//sum
				++peak_count;
				intensity_sum += i->getIntensity();
  			mzint[int(i->getMZ()*mult)] += i->getIntensity();
				mzcount[int(i->getMZ()*mult)]++;
				mzsum[int(i->getMZ()*mult)] += i->getMZ();

				rt[i.getRT()] += i->getIntensity();
				//max
				intensity_max = max(intensity_max,(DoubleReal)(i->getIntensity()));
			}
		}


		// write to spectra
		MSExperiment<>::SpectrumType::ContainerType& cont_mz = projection_mz_[0].getContainer();
		MSExperiment<>::SpectrumType::ContainerType& cont_rt = projection_rt_[0].getContainer();
		
		//resize and add boundary peaks		
		cont_mz.resize(mzint.size()+2);
		cont_mz[0].setMZ(visible_area_.min()[0]);
		cont_mz[0].setIntensity(0.0);
		cont_mz[1].setMZ(visible_area_.max()[0]);
		cont_mz[1].setIntensity(0.0);
		cont_rt.resize(rt.size()+2);
		cont_rt[0].setMZ(visible_area_.min()[1]);
		cont_rt[0].setIntensity(0.0);
		cont_rt[1].setMZ(visible_area_.max()[1]);
		cont_rt[1].setIntensity(0.0);
		
		UInt i = 2;
		map<int,float>::iterator intit = mzint.begin();
		map<int,int>::iterator cit = mzcount.begin();
		
		for (map<int, float>::iterator it = mzsum.begin(); it != mzsum.end(); ++it)
		{
			cont_mz[i].setMZ(it->second/cit->second);
			cont_mz[i].setIntensity(intit->second);
			intit++;
			cit++;
			++i;
		}

		i = 2;
		for (map<float, float>::iterator it = rt.begin(); it != rt.end(); ++it)
		{
			cont_rt[i].setMZ(it->first);
			cont_rt[i].setIntensity(it->second);
			++i;
		}
		
		if (isMzToXAxis())
		{
			emit showProjectionHorizontal(projection_mz_,Spectrum1DCanvas::DM_PEAKS);
			emit showProjectionVertical(projection_rt_,Spectrum1DCanvas::DM_CONNECTEDLINES);	
		}
		else
		{
			emit showProjectionHorizontal(projection_rt_,Spectrum1DCanvas::DM_CONNECTEDLINES);
			emit showProjectionVertical(projection_mz_,Spectrum1DCanvas::DM_PEAKS);
		}
		showProjectionInfo(peak_count, intensity_sum, intensity_max);
	}
	
	Int Spectrum2DCanvas::finishAdding(float low_intensity_cutoff)
	{
		current_layer_ = getLayerCount()-1;
		
		recalculateDotGradient_(activeLayerIndex());
		
		if (layers_.back().type==LayerData::DT_PEAK) //peak data
		{
			currentPeakData_().sortSpectra(true);
			currentPeakData_().updateRanges(1);
			getCurrentLayer_().min_int = low_intensity_cutoff;
			getCurrentLayer_().max_int = getCurrentPeakData().getMaxInt();
			update_buffer_ = true;
		}
		else //feature data
		{
			getCurrentLayer_().features.updateRanges();
			getCurrentLayer_().max_int = getCurrentLayer().features.getMaxInt();
			setLayerFlag(LayerData::F_HULLS,true);
		}
		
		//overall values update
		recalculateRanges_(0,1,2);
		//cout << "New data range: " << overall_data_range_ << endl;
		
		if (getLayerCount()==1)
		{
			AreaType tmp_area;
			tmp_area.assign(overall_data_range_);
			visible_area_ = tmp_area;
			emit visibleAreaChanged(tmp_area);
		}
		else
		{
			resetZoom();
		}
		
		if (getLayerCount()==2)
		{
			setIntensityMode(IM_PERCENTAGE);
		}
		
		intensityModeChange_();
		
		emit sendStatusMessage("",0);
		
		emit layerActivated(this);

		return current_layer_;
	}
	
	void Spectrum2DCanvas::removeLayer(int layer_index )
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()))
		{
			return;
		}

		//unselect all peaks
		delete(selected_peak_);
		selected_peak_ = 0;
		delete(measurement_start_);
		measurement_start_ = 0;
		delete(measurement_stop_);
		measurement_stop_ = 0;
	
		//remove the data
		layers_.erase(layers_.begin()+layer_index);
		
		//update visible area and boundaries
		recalculateRanges_(0,1,2);

		AreaType tmp;
		tmp.assign(overall_data_range_);
		if (tmp != visible_area_)
		{
			visible_area_.assign(overall_data_range_);
		}

		//update current layer
		if (current_layer_!=0 && current_layer_ >= getLayerCount())
		{
			current_layer_ = getLayerCount()-1;
		}

		if (layers_.empty())
		{
			return;
		}

		intensityModeChange_();
		emit layerActivated(this);
	}
	
	//change the current layer
	void Spectrum2DCanvas::activateLayer(int layer_index)
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()) || layer_index==int(current_layer_))
		{
			return ;
		}
		current_layer_ = layer_index;
		emit layerActivated(this);
		
		//unselect all peaks
		delete(selected_peak_);
		selected_peak_ = 0;
		delete(measurement_start_);
		measurement_start_ = 0;
		delete(measurement_stop_);
		measurement_stop_ = 0;
		
		update();
	}

	void Spectrum2DCanvas::recalculateSnapFactor_()
	{
		if (intensity_mode_ == IM_SNAP) 
		{
			double local_max  = -numeric_limits<double>::max();
			for (UInt i=0; i<getLayerCount(); i++)
			{
				if (getLayer(i).visible)
				{
					if (getLayer(i).type==LayerData::DT_PEAK)
					{
						for (ExperimentType::ConstAreaIterator it = getPeakData(i).areaBeginConst(visible_area_.min()[1],visible_area_.max()[1],visible_area_.min()[0],visible_area_.max()[0]); 
								 it != getPeakData(i).areaEndConst(); 
								 ++it)
						{
							if (it->getIntensity() > local_max && it->getIntensity()>=getLayer(i).min_int && it->getIntensity()<=getLayer(i).max_int)
							{
								local_max = it->getIntensity();
							}
						}
					}
					else //features
					{
						for (FeatureMapType::ConstIterator it = getLayer(i).features.begin();
							   it != getLayer(i).features.end();
							   ++it)
						{
							if ( it->getRT() >= visible_area_.min()[1] &&
									 it->getRT() <= visible_area_.max()[1] &&
									 it->getMZ() >= visible_area_.min()[0] &&
									 it->getMZ() <= visible_area_.max()[0] &&
									 it->getIntensity()>=getLayer(i).min_int && 
									 it->getIntensity()<=getLayer(i).max_int &&
									 it->getIntensity() > local_max)
							{
								local_max = it->getIntensity();
							}
						}
					}
				}
			}
			snap_factor_ = overall_data_range_.max()[2]/local_max;			
		}
		else
		{ 
			snap_factor_ = 1.0;
		}
	}

	void Spectrum2DCanvas::updateScrollbars_()
	{
		if (isMzToXAxis())
		{
			emit updateHScrollbar(overall_data_range_.min()[0],visible_area_.min()[0],visible_area_.max()[0],overall_data_range_.max()[0]);
			emit updateVScrollbar(overall_data_range_.min()[1],visible_area_.min()[1],visible_area_.max()[1],overall_data_range_.max()[1]);
		}
		else
		{
			
			emit updateVScrollbar(overall_data_range_.min()[0],visible_area_.min()[0],visible_area_.max()[0],overall_data_range_.max()[0]);
			emit updateHScrollbar(overall_data_range_.min()[1],visible_area_.min()[1],visible_area_.max()[1],overall_data_range_.max()[1]);
		}
	}

	void Spectrum2DCanvas::horizontalScrollBarChange(int value)
	{
		AreaType new_area = visible_area_;
		if (isMzToXAxis())
		{
			new_area.setMinX(value);
			new_area.setMaxX(value + (visible_area_.maxX() - visible_area_.minX()));
			//cout << __PRETTY_FUNCTION__ << endl;
			changeVisibleArea_(new_area);
		}
		else
		{
			new_area.setMinY(value);
			new_area.setMaxY(value + (visible_area_.maxY() - visible_area_.minY()));
			//cout << __PRETTY_FUNCTION__ << endl;
			changeVisibleArea_(new_area);
		}
	}

	void Spectrum2DCanvas::verticalScrollBarChange(int value)
	{
		AreaType new_area = visible_area_;
		if (!isMzToXAxis())
		{
			new_area.setMinX(value);
			new_area.setMaxX(value + (visible_area_.maxX() - visible_area_.minX()));
			//cout << __PRETTY_FUNCTION__ << endl;
			changeVisibleArea_(new_area);
		}
		else
		{
			new_area.setMinY(value);
			new_area.setMaxY(value + (visible_area_.maxY() - visible_area_.minY()));
			//cout << __PRETTY_FUNCTION__ << endl;
			changeVisibleArea_(new_area);
		}
	}

	void Spectrum2DCanvas::paintEvent(QPaintEvent* e)
	{
#ifdef DEBUG_TOPPVIEW
		cout << "BEGIN " << __PRETTY_FUNCTION__ << endl;
	  cout << "  Visible area -- m/z: " << visible_area_.minX() << " - " << visible_area_.maxX() << " int: " << visible_area_.minY() << " - " << visible_area_.maxY() << endl;
	  cout << "  Overall area -- m/z: " << overall_data_range_.min()[0] << " - " << overall_data_range_.max()[0] << " int: " << overall_data_range_.min()[1] << " - " << overall_data_range_.max()[1] << endl; 
#endif
#ifdef TIMING_TOPPVIEW
		QTime timer;
 		timer.start();
#endif
		
		QPainter painter;
		
		if (update_buffer_)
		{
			update_buffer_ = false;
			
			//recalculate snap factor
			recalculateSnapFactor_();
			
			buffer_.fill(QColor(param_.getValue("background_color").toQString()).rgb());
			painter.begin(&buffer_);

			for (UInt i=0; i<getLayerCount(); i++)
			{
				if (getLayer(i).visible)
				{
					if (getLayer(i).type==LayerData::DT_PEAK)
					{
						paintDots_(i, painter);
					}
					else if (getLayer(i).type==LayerData::DT_FEATURE)
					{
						//cout << "dot feature layer: " << i << endl;
						paintDots_(i, painter);
						if (getLayerFlag(i,LayerData::F_HULLS))
						{
							paintConvexHulls_(i, painter);
						}
					}
					else if (getLayer(i).type==LayerData::DT_FEATURE_PAIR)
					{
						//cout << "dot feature pair layer: " << i << endl;
						paintDots_(i, painter);
						if( getLayerFlag(i,LayerData::F_HULLS))
						{
							paintConvexHulls_(i, painter);
						}
						paintFeaturePairConnections_(i, painter);
					}
				}
			}
			paintGridLines_(painter);
			painter.end();
		}
		
		painter.begin(this);
		
		//copy peak data from buffer
		QVector<QRect> rects = e->region().rects();
		for (int i = 0; i < (int)rects.size(); ++i)
		{
			painter.drawPixmap(rects[i].topLeft(), buffer_, rects[i]);
		}
		
		//draw mesaurement peak
		if (measurement_start_)
		{
			painter.setPen(Qt::black);
			
			QPoint line_begin, line_end;
			
			if (measurement_stop_)
			{
				 dataToWidget_(measurement_stop_->getMZ(), measurement_stop_->getRT(), line_end);
				//cout << "Line end: " << line_end << endl;
			}
			else
			{
				line_end = last_mouse_pos_;
				//cout << "Ende: " << line_end.x() << " " << line_end.y() << endl;
			}
			dataToWidget_(measurement_start_->getMZ(), measurement_start_->getRT(), line_begin);
			painter.drawLine(line_begin, line_end);
		}
		highlightPeak_(painter, measurement_start_);
		highlightPeak_(painter, measurement_stop_);
		
		//draw selected peak
		highlightPeak_(painter, selected_peak_);
		
		//draw convex hull of selected peak
		if (selected_peak_ && getCurrentLayer().type!=LayerData::DT_PEAK)
		{
			painter.setPen(QPen(Qt::red, 2));
			paintConvexHulls_(selected_peak_->getConvexHulls(),painter);
		}
		
		painter.end();
#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
#ifdef TIMING_TOPPVIEW
		cout << "2D PaintEvent took " << timer.elapsed() << " ms" << endl << endl;
#endif	
	}

	void Spectrum2DCanvas::mousePressEvent(QMouseEvent* e)
	{
		last_mouse_pos_ = e->pos();
		
		switch (action_mode_)
		{
			case AM_SELECT:
				if (e->button() == Qt::LeftButton)
				{
					if (e->modifiers() & Qt::ControlModifier) //measure
					{
						if (selected_peak_)
						{
							delete(measurement_start_);
							measurement_start_ = new Feature(*selected_peak_);
						}
						else
						{
							delete(measurement_start_);
							measurement_start_ = 0;
						}
						delete(measurement_stop_);
						measurement_stop_ = 0;
					}
				}
				break;

			case AM_ZOOM:
				if (e->button() == Qt::LeftButton)
				{
					if (e->modifiers() & Qt::ControlModifier) //translate
					{
						setCursor(Qt::ClosedHandCursor);
					}
					else //zoom
					{
						rubber_band_.setGeometry(e->pos().x(),e->pos().y(),0,0);
						rubber_band_.show();
					}
				}
				break;
			default:
				break;
		}

		e->accept();
	}

	void Spectrum2DCanvas::mouseMoveEvent(QMouseEvent* e)
	{
		QPoint pos = e->pos();
		
		switch (action_mode_)
		{
			case AM_SELECT:
			{
				if (e->modifiers() & Qt::ControlModifier) //measure
				{
					// highlight nearest peak
					if (e->buttons() == Qt::NoButton)
					{
						Feature* max_peak = findNearestPeak_(pos);
						
						if (max_peak && max_peak != selected_peak_ && !measurement_start_)
						{
							//show Peak Coordinates
							emit sendCursorStatus(max_peak->getMZ(), max_peak->getIntensity(), max_peak->getRT());
							//show status message (label + charge)
							String status;
							String label = max_peak->getMetaValue(3).toString();
							if (label!="") status = status + " Label: " + label;
							String charge = max_peak->getCharge();
							if (charge!="") status = status + " Charge: " + charge;
							if (status!="") sendStatusMessage(status, 0);
						}
						
						selected_peak_ = max_peak;
						update();
					}
					else if (e->buttons() & Qt::LeftButton && measurement_start_)
					{
						measurement_stop_ = findNearestPeak_(pos);
						last_mouse_pos_ = pos;
						update();
					
						if (measurement_stop_)
						{
							emit sendCursorStatus(measurement_stop_->getMZ() - measurement_start_->getMZ(),
									      measurement_stop_->getIntensity() / measurement_start_->getIntensity(),
									      measurement_stop_->getRT() - measurement_start_->getRT());
						}
						else
						{
							emit sendCursorStatus(measurement_start_->getMZ(), measurement_start_->getIntensity(), measurement_start_->getRT());
						}
					}
				}
				else //select 	 
				{ 	 
					// highlight nearest peak 	 
					if (e->buttons() == Qt::NoButton) 	 
					{ 	 
						Feature* max_peak = findNearestPeak_(pos); 	 
						if (max_peak) 	 
						{ 	 
							// show Peak Coordinates (with intensity) 	 
							emit sendCursorStatus(max_peak->getMZ(), max_peak->getIntensity(), max_peak->getRT()); 	  
							//show status message (label + charge)
							String status;
							String label = max_peak->getMetaValue(3).toString();
							if (label!="") status = status + " Label: " + label;
							String charge = max_peak->getCharge();
							if (charge!="") status = status + " Charge: " + charge;
							if (status!="") sendStatusMessage(status, 0);	 
						} 	 
						else 	 
						{ 	 
						 //show Peak Coordinates (without intensity) 	 
						 PointType pnt = widgetToData_(pos); 	 
						 emit sendCursorStatus( pnt[0], -1.0, pnt[1]); 	 
						} 	 
						
						selected_peak_ = max_peak; 	 
						update(); 	 
					}
				}
				break;
			}
			case AM_ZOOM:
			{
				//show Peak Coordinates
				PointType pnt = widgetToData_(pos);
				emit sendCursorStatus( pnt[0], -1.0, pnt[1]);
				
				if (e->buttons() & Qt::LeftButton)
				{
					if (e->modifiers() & Qt::ControlModifier) //translate
					{
						//caldulate data coordinates of shift
						PointType old_data = widgetToData_(last_mouse_pos_);
						PointType new_data = widgetToData_(pos);
						//calculate x shift
						double shift = old_data.getX() - new_data.getX();
						double newLoX = visible_area_.minX() + shift;
						double newHiX = visible_area_.maxX() + shift;
						// check if we are falling out of bounds
						if (newLoX < overall_data_range_.minX())
						{
							newLoX = overall_data_range_.minX();
							newHiX = newLoX + visible_area_.width();
						}
						if (newHiX > overall_data_range_.maxX())
						{
							newHiX = overall_data_range_.maxX();
							newLoX = newHiX - visible_area_.width();
						}
						//calculate y shift
						shift = old_data.getY() - new_data.getY();
						double newLoY = visible_area_.minY() + shift;
						double newHiY = visible_area_.maxY() + shift;
						// check if we are falling out of bounds
						if (newLoY < overall_data_range_.minY())
						{
							newLoY = overall_data_range_.minY();
							newHiY = newLoY + visible_area_.height();
						}
						if (newHiY > overall_data_range_.maxY())
						{
							newHiY = overall_data_range_.maxY();
							newLoY = newHiY - visible_area_.height();
						}
		     		
		     		//change area
						//cout << "New area: x " << newLoX <<"-"<< newHiX << " - y "<<newLoY <<"-"<< newHiY << endl;
						//cout << __PRETTY_FUNCTION__ << endl;
						changeVisibleArea_(AreaType(newLoX,newLoY,newHiX,newHiY));
		
						last_mouse_pos_ = pos;
					}
					else //zoom
					{
						rubber_band_.setGeometry(last_mouse_pos_.x(), last_mouse_pos_.y(), pos.x() - last_mouse_pos_.x(), pos.y() - last_mouse_pos_.y());
						update();
					}
				}
				else
				{
					if (e->modifiers() & Qt::ControlModifier) //translate
					{
						setCursor(Qt::OpenHandCursor);
					}
					else //zoom
					{
						setCursor(Qt::CrossCursor);
					}
				}
				break;
			}
			default:
				break;
		}
		e->accept();
	}
	
	void Spectrum2DCanvas::mouseReleaseEvent(QMouseEvent* e)
	{
		QPoint pos = e->pos();

		switch (action_mode_)
		{
			case AM_SELECT:
			{
				if (e->button() == Qt::LeftButton)
				{
					if (e->modifiers() & Qt::ControlModifier) //measure
					{
						if (!measurement_stop_)
						{
							delete(measurement_start_);
							measurement_start_ = 0;
						}
						else
						{
							measurement_stop_ = new Feature(*measurement_stop_);
						}
						
						update();
						
						if (measurement_start_)
						{
							emit sendStatusMessage(QString("Measured: dRT = %1, dMZ = %3, Intensity ratio = %2")
																		.arg(measurement_stop_->getRT() - measurement_start_->getRT())
																		.arg(measurement_stop_->getIntensity() / measurement_start_->getIntensity())
																		.arg(measurement_stop_->getMZ() - measurement_start_->getMZ()).toAscii().data(), 0);
						}

					}
				}
				break;
			}
			case AM_ZOOM:
			{
				if(e->button() == Qt::LeftButton)
				{
					rubber_band_.hide();
					if (e->modifiers() & Qt::ControlModifier) //translate
					{
						setCursor(Qt::OpenHandCursor);
					}
					else //zoom
					{
						QRect rect = rubber_band_.geometry();
						if (rect.width()!=0 && rect.height()!=0) //probably double-click -> mouseDoubleClickEvent
						{
							AreaType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
							//cout << __PRETTY_FUNCTION__ << endl;
							changeVisibleArea_(area, true);
						}
					}
				}
				break;
			}
			default:
				break;
		}
		e->accept();
	}

	void Spectrum2DCanvas::mouseDoubleClickEvent(QMouseEvent* e)
	{
		// left-doubleclick shows the whole spectrum
		if (e->button() == Qt::LeftButton && action_mode_ == AM_ZOOM)
		{
			resetZoom();
		}
	}

	void Spectrum2DCanvas::wheelEvent(QWheelEvent* e)
	{
		if (e->delta() > 0) // forward rotation -> zoom in
		{
			PointType new_pos = visible_area_.center();
			float half_width = visible_area_.width() / 2.0 * 0.9;
			float half_height = visible_area_.height() / 2.0f * 0.9;
			
			//cout << __PRETTY_FUNCTION__ << endl;
			changeVisibleArea_(AreaType(new_pos.getX() - half_width, new_pos.getY() - half_height, new_pos.getX() + half_width, new_pos.getY() + half_height), true);
		}
		else // backward rotation -> zoom out
		{
			zoomBack_();
		}
		e->accept();
	}

	void Spectrum2DCanvas::contextMenuEvent(QContextMenuEvent* e)
	{
		DoubleReal rt = widgetToData_(e->pos())[1];
		
		const LayerData& layer = getCurrentLayer();

		QMenu context_menu(this);
		QAction* a = 0;
		QAction* result = 0;
		
		//-------------------PEAKS----------------------------------
		if (layer.type==LayerData::DT_PEAK)
		{
			//select surrounding scans
			MSExperiment<>::ConstIterator first = getCurrentPeakData().begin();
			MSExperiment<>::ConstIterator it = getCurrentPeakData().RTBegin(rt);
			MSExperiment<>::ConstIterator begin = it;
			MSExperiment<>::ConstIterator end = it;
			UInt count = 0;
			while (begin!=first && count <4)
			{
				--begin;
				++count;
			}
			count = 0;
			while (end!=getCurrentPeakData().end() && count <5)
			{
				++end;
				++count;
			}
			
			//add scans
			QMenu* scans = context_menu.addMenu("View MS scan in 1D");
			QMenu* meta = context_menu.addMenu("View/edit meta data");
			while(begin!=end)
			{
				if (begin==it)
				{
					scans->addSeparator();
					meta->addSeparator();
				}
				if (begin->getMSLevel()<2)
				{
					a = scans->addAction(QString("RT: ") + QString::number(begin->getRT()));
					a->setData((int)(begin-first));
					a = meta->addAction(QString("RT: ") + QString::number(begin->getRT()));
					a->setData((int)(begin-first));
				}
				else
				{
					a = scans->addAction(QString("RT: ") + QString::number(begin->getRT()) + "  Precursor m/z:" + QString::number(begin->getPrecursorPeak().getPosition()[0]));
					a->setData((int)(begin-first));
					a = meta->addAction(QString("RT: ") + QString::number(begin->getRT()) + "  Precursor m/z:" + QString::number(begin->getPrecursorPeak().getPosition()[0]));
					a->setData((int)(begin-first));
				}
				if (begin==it)
				{
					scans->addSeparator();
					meta->addSeparator();
				}
				++begin;
			}
			context_menu.addAction("View visible data in 3D");
			
			if ((result = context_menu.exec(mapToGlobal(e->pos()))))
			{
				if (result->parent()==scans)
				{
					emit showSpectrumAs1D(result->data().toInt());
				}
				else if (result->parent()==meta)
				{
					MSMetaDataExplorer dlg(true, this);
		      dlg.setWindowTitle("View/Edit meta data");
		    	dlg.visualize(static_cast<SpectrumSettings&>(currentPeakData_()[result->data().toInt()]));
		      dlg.exec();
				}
				else if (result->text() == "View visible data in 3D")
				{
					emit showCurrentPeaksAs3D();
				}
			}			
		}
		//-------------------FEATURES----------------------------------
		else if (layer.type==LayerData::DT_FEATURE || layer.type==LayerData::DT_FEATURE_PAIR)
		{
			//search for nearby features
			DPosition<2> p1 = widgetToData_(e->pos()+ QPoint(10,10));
			DPosition<2> p2   = widgetToData_(e->pos()- QPoint(10,10));
			DoubleReal rt_min = min(p1[1],p2[1]);
			DoubleReal rt_max = max(p1[1],p2[1]);
			DoubleReal mz_min = min(p1[0],p2[0]);
			DoubleReal mz_max = max(p1[0],p2[0]);
			
			QMenu* meta = context_menu.addMenu("View/edit meta data");
			bool present = false;
			FeatureMapType& features = getCurrentLayer_().features;
			for (FeatureMapType::Iterator it = features.begin(); it!=features.end(); ++it)
			{
				if (it->getMZ() <= mz_max && it->getMZ() >= mz_min && it->getRT() <= rt_max && it->getRT() >= rt_min)
				{
					present=true;
					a = meta->addAction(QString("RT: ") + QString::number(it->getRT()) + "  m/z:" + QString::number(it->getMZ()));
					a->setData((int)(it-features.begin()));
				}  
			}
			
			if (present && (result = context_menu.exec(mapToGlobal(e->pos()))))
			{
				MSMetaDataExplorer dlg(true, this);
	      dlg.setWindowTitle("View/Edit meta data");
	      
	      vector<PeptideIdentification>& ids = features[result->data().toInt()].getPeptideIdentifications();
				for (vector<PeptideIdentification>::iterator it=ids.begin(); it!=ids.end(); ++it)
				{
	    		dlg.visualize(*it);
				}
				
	      dlg.exec();
			}
		}
		
		e->accept();
	}

	void Spectrum2DCanvas::showCurrentLayerPreferences()
	{
		Internal::Spectrum2DPrefDialog dlg(this);

		ColorSelector* bg_color = dlg.findChild<ColorSelector*>("bg_color");
		QComboBox* mapping = dlg.findChild<QComboBox*>("mapping");
		MultiGradientSelector* gradient = dlg.findChild<MultiGradientSelector*>("gradient");

		bg_color->setColor(QColor(param_.getValue("background_color").toQString()));
		if (isMzToXAxis())
		{
			mapping->setCurrentIndex(0);
		}
		else
		{
			mapping->setCurrentIndex(1);
		}
		gradient->gradient().fromString(getCurrentLayer_().param.getValue("dot:gradient"));
		
		if (dlg.exec())
		{
			param_.setValue("background_color",bg_color->getColor().name().toAscii().data());
			if ((mapping->currentIndex()==0 && !isMzToXAxis()) || (mapping->currentIndex()==1 && isMzToXAxis()))
			{
				mzToXAxis(!isMzToXAxis());
			}
			getCurrentLayer_().param.setValue("dot:gradient",gradient->gradient().toString());
			
			recalculateDotGradient_(activeLayerIndex());
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
	}

} //namespace OpenMS

