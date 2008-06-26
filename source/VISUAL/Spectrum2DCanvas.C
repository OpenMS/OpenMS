// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DPrefDialog.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/SYSTEM/FileWatcher.h>

//STL
#include <algorithm>	

//QT
#include <QtGui/QMouseEvent>
#include <QtGui/QPainter>
#include <QtGui/QMenu>
#include <QtGui/QBitmap>
#include <QtGui/QPolygon>
#include <QtCore/QTime>
#include <QtGui/QComboBox>
#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>

using namespace std;

namespace OpenMS
{
	using namespace Internal;

	Spectrum2DCanvas::Spectrum2DCanvas(const Param& preferences, QWidget* parent)
		: SpectrumCanvas(preferences, parent),
			selected_peak_(),
			measurement_start_()
	{
    //Paramater handling
    defaults_.setValue("background_color", "#ffffff", "Background color.");
    defaults_.setValue("interpolation_steps", 200, "Number of interploation steps for peak gradient precalculation.");
    defaults_.setMinInt("interpolation_steps",1);
    defaults_.setMaxInt("interpolation_steps",1000);
    defaults_.setValue("dot:gradient", "Linear|0,#ffea00;6,#ff0000;14,#aa00ff;23,#5500ff;100,#000000", "Multi-color gradient for peaks.");
    defaults_.setValue("mapping_of_mz_to","x_axis","Determines with axis is the m/z axis.");
		vector<String> strings;
		strings.push_back("x_axis");
		strings.push_back("y_axis");
		defaults_.setValidStrings("mapping_of_mz_to",strings);
		defaultsToParam_();
		setName("Spectrum2DCanvas");
		setParameters(preferences);

		projection_mz_.resize(1);
		projection_rt_.resize(1);
		
		//set preferences and update widgets acoordningly
		if (String(param_.getValue("mapping_of_mz_to")) != "x_axis")
		{
			mzToXAxis(false);
		}
	}
	
	Spectrum2DCanvas::~Spectrum2DCanvas()
	{
		//cout << "DEST Spectrum2DCanvas" << endl;
	}
	
	void Spectrum2DCanvas::highlightPeak_(QPainter& painter, const PeakIndex& peak)
	{
		if (!peak.isValid()) return;
		painter.save();
		painter.setPen(QPen(Qt::red, 2));
		QPoint pos;
		if (getCurrentLayer().type!=LayerData::DT_PEAK)
		{
			dataToWidget_(peak.getFeature(getCurrentLayer().features).getMZ(),peak.getFeature(getCurrentLayer().features).getRT(),pos);
		}
		else
		{
			dataToWidget_(peak.getPeak(getCurrentLayer().peaks).getMZ(),peak.getSpectrum(getCurrentLayer().peaks).getRT(),pos);
		}
		painter.drawEllipse(pos.x() - 5, pos.y() - 5, 10, 10);
		painter.restore();
	}
	
	PeakIndex Spectrum2DCanvas::findNearestPeak_(const QPoint& pos)
	{
		///no layers => return invalid peak index
		if (layers_.empty()) return PeakIndex();
		
		//Constructing the area corrects swapped mapping of RT and m/z
		AreaType area (widgetToData_(pos - QPoint(5,5)),widgetToData_(pos + QPoint(5,5)));

		float max_int = -1 * numeric_limits<float>::max();

		if (getCurrentLayer().type==LayerData::DT_PEAK)
		{
			for (ExperimentType::ConstAreaIterator i = getCurrentLayer().peaks.areaBeginConst(area.min()[1],area.max()[1],area.min()[0],area.max()[0]); 
					 i != getCurrentLayer().peaks.areaEndConst(); 
					 ++i)
			{
				PeakIndex pi = i.getPeakIndex();
				if (i->getIntensity() > max_int && getCurrentLayer().filters.passes(getCurrentLayer().peaks[pi.spectrum],pi.peak))
				{
					//cout << "new max: " << i.getRT() << " " << i->getMZ() << endl;
					max_int = i->getIntensity();
					return i.getPeakIndex();
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
						 getCurrentLayer().filters.passes(*i) )
				{
					if (i->getIntensity() > max_int)
					{
						max_int = i->getIntensity();
												
						return PeakIndex(i-getCurrentLayer().features.begin());
					}				
				}
			}	 	
		}
		return PeakIndex();
	}
	
	float Spectrum2DCanvas::betweenFactor_(float v1, float v2, float val)
	{
		float lo = min(v1, v2);
		float hi = max(v1, v2);
		return (hi - lo == 0) ? 1 : (val - lo) / (hi - lo);
	}

	void Spectrum2DCanvas::paintDots_(UInt layer_index, QPainter& painter)
	{
#ifdef TIMING_TOPPVIEW
		QTime timer;
		timer.start();
#endif	
		
		const LayerData& layer = getLayer(layer_index);
		
		percentage_factor_ = 1.0;
		if (intensity_mode_ == IM_PERCENTAGE)
		{
			if (layer.type == LayerData::DT_PEAK && layer.peaks.getMaxInt()>0.0)
			{
				percentage_factor_ = overall_data_range_.max()[2]/layer.peaks.getMaxInt();
			}
			else if (/* layer.type != LayerData::DT_PEAK && */ layer.features.getMaxInt()>0.0)
			{
				percentage_factor_ = overall_data_range_.max()[2]/layer.features.getMaxInt();
			}
		}
		
		//temporary variable
		QPoint pos;
		
		painter.setPen(Qt::black);
		
		if (layer.type==LayerData::DT_PEAK) //peaks
		{
			//determine number of MS1 scans
			UInt scans = 0;
			ExperimentType::ConstIterator it = layer.peaks.RTBegin(visible_area_.min()[1]);
			while (it != layer.peaks.RTEnd(visible_area_.max()[1]))
			{
				if (it->getMSLevel()==1) ++scans;
				++it;
			}
			//determine number of shown peaks
			Int peaks = 0;
			it = layer.peaks.RTBegin(visible_area_.min()[1]) + scans/2;
			while (it!=layer.peaks.end() && it->getMSLevel()!=1)
			{
				++it;
			}
			if (it!=layer.peaks.end())
			{
				ExperimentType::SpectrumType::ConstIterator it2 = it->MZBegin(visible_area_.min()[0]);
				while (it2!=it->MZEnd(visible_area_.max()[0]))
				{
					if (layer.filters.passes(*it,it2-it->begin())) ++peaks;
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
			for (ExperimentType::ConstAreaIterator i = layer.peaks.areaBeginConst(visible_area_.min()[1],visible_area_.max()[1],visible_area_.min()[0],visible_area_.max()[0]); 
					 i != layer.peaks.areaEndConst(); 
					 ++i)
			{
				PeakIndex pi = i.getPeakIndex();
				if (layer.filters.passes(layer.peaks[pi.spectrum],pi.peak))
				{
					painter.setPen(heightColor_(i->getIntensity(), layer.gradient));
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
				const ExperimentType& exp = layer.peaks; 
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
			for (FeatureMapType::ConstIterator i = layer.features.begin();
				   i != layer.features.end();
				   ++i)
			{
				if ( i->getRT() >= visible_area_.min()[1] &&
						 i->getRT() <= visible_area_.max()[1] &&
						 i->getMZ() >= visible_area_.min()[0] &&
						 i->getMZ() <= visible_area_.max()[0] &&
						 layer.filters.passes(*i))
				{
					painter.setPen(heightColor_(i->getIntensity(), layer.gradient));
					dataToWidget_(i->getMZ(),i->getRT(),pos);
					painter.drawLine(pos.x(),pos.y()-1,pos.x(),pos.y()+1);
					painter.drawLine(pos.x()-1,pos.y(),pos.x()+1,pos.y());
					if (numbers)
					{
						painter.setPen(Qt::black);
						//paint label of feature number
						QString label = QString::number(num);
						if (i->metaValueExists(3))
						{
							label.append(" (").append(i->getMetaValue(3).toString().c_str()).append(")");
						}
						painter.drawText(pos.x()+10,pos.y()+10,label);
					}
				}
				++num;
			}
		}

#ifdef TIMING_TOPPVIEW
		cout << "paintDots_ took " << timer.elapsed() << " ms" << endl;
#endif	
	}

	void Spectrum2DCanvas::paintTraceConvexHulls_(UInt layer_index, QPainter& painter)
	{
		painter.setPen(Qt::black);
		
		const LayerData& layer = getLayer(layer_index);
		for (FeatureMapType::ConstIterator i = layer.features.begin(); i != layer.features.end(); ++i)
		{
			if ( i->getRT() >= visible_area_.min()[1] &&
					 i->getRT() <= visible_area_.max()[1] &&
					 i->getMZ() >= visible_area_.min()[0] &&
					 i->getMZ() <= visible_area_.max()[0] &&
					 layer.filters.passes(*i)
				 )
			{
				paintConvexHulls_(i->getConvexHulls(),painter);
			}
		}
	}            

	void Spectrum2DCanvas::paintFeatureConvexHulls_(UInt layer_index, QPainter& painter)
	{
		painter.setPen(Qt::black);
		const LayerData& layer = getLayer(layer_index);
		for (FeatureMapType::ConstIterator i = layer.features.begin(); i != layer.features.end(); ++i)
		{
			if ( i->getRT() >= visible_area_.min()[1] &&
					 i->getRT() <= visible_area_.max()[1] &&
					 i->getMZ() >= visible_area_.min()[0] &&
					 i->getMZ() <= visible_area_.max()[0] &&
					 layer.filters.passes(*i))
			{				
				//paint hull points
				ConvexHull2D hull = i->getConvexHull();
				QPolygon points;
				points.resize(hull.getPoints().size());
		
				UInt index=0;
				QPoint pos;
				//iterate over hull points
				for(ConvexHull2D::PointArrayType::const_iterator it=hull.getPoints().begin(); it!=hull.getPoints().end(); ++it, ++index)
				{
					dataToWidget_(it->getY(), it->getX(),pos);
					points.setPoint(index, pos);
				}	
				//cout << "Hull: " << hull << " Points: " << points.size()<<endl;
				painter.drawPolygon(points);
			}
		}
	}    

  void Spectrum2DCanvas::paintConvexHulls_(const vector<ConvexHull2D>& hulls, QPainter& painter)
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
		getLayer_(layer).gradient.activatePrecalculationMode(min(0.0,overall_data_range_.min()[2]), overall_data_range_.max()[2], param_.getValue("interpolation_steps"));
	}
	
	void Spectrum2DCanvas::updateProjections()
	{
		const LayerData& layer = getCurrentLayer();
		if (layer.type != LayerData::DT_PEAK)
		{
			QMessageBox::critical(this,"Error","Cannot show projections of feature layers!");
			return;	
		}
		
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


		for (ExperimentType::ConstAreaIterator i = layer.peaks.areaBeginConst(visible_area_.min()[1],visible_area_.max()[1],visible_area_.min()[0],visible_area_.max()[0]); 
				 i != layer.peaks.areaEndConst();
				 ++i)
		{
			PeakIndex pi = i.getPeakIndex();
			if (layer.filters.passes(layer.peaks[pi.spectrum],pi.peak))
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
		projection_mz_[0].resize(mzint.size()+2);
		projection_mz_[0][0].setMZ(visible_area_.min()[0]);
		projection_mz_[0][0].setIntensity(0.0);
		projection_mz_[0][1].setMZ(visible_area_.max()[0]);
		projection_mz_[0][1].setIntensity(0.0);
		projection_rt_[0].resize(rt.size()+2);
		projection_rt_[0][0].setMZ(visible_area_.min()[1]);
		projection_rt_[0][0].setIntensity(0.0);
		projection_rt_[0][1].setMZ(visible_area_.max()[1]);
		projection_rt_[0][1].setIntensity(0.0);
		
		UInt i = 2;
		map<int,float>::iterator intit = mzint.begin();
		map<int,int>::iterator cit = mzcount.begin();
		
		for (map<int, float>::iterator it = mzsum.begin(); it != mzsum.end(); ++it)
		{
			projection_mz_[0][i].setMZ(it->second/cit->second);
			projection_mz_[0][i].setIntensity(intit->second);
			intit++;
			cit++;
			++i;
		}

		i = 2;
		for (map<float, float>::iterator it = rt.begin(); it != rt.end(); ++it)
		{
			projection_rt_[0][i].setMZ(it->first);
			projection_rt_[0][i].setIntensity(it->second);
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
	
	bool Spectrum2DCanvas::finishAdding_()
	{
		//unselect all peaks
		selected_peak_.clear();
		measurement_start_.clear();
		
		current_layer_ = getLayerCount()-1;
		
		if (layers_.back().type==LayerData::DT_PEAK) //peak data
		{
			currentPeakData_().sortSpectra(true);
			currentPeakData_().updateRanges(1);
			
			update_buffer_ = true;
			
			//Abort if no data points are contained
			if (currentPeakData_().size()==0 || currentPeakData_().getSize()==0)
			{
				layers_.resize(getLayerCount()-1);
				if (current_layer_!=0) current_layer_ = current_layer_-1;
				QMessageBox::critical(this,"Error","Cannot add a dataset that contains no survey scans. Aborting!");
				return false;
			}
		}
		else //feature data
		{
			getCurrentLayer_().features.updateRanges();
			setLayerFlag(LayerData::F_HULL,true);

			//Abort if no data points are contained
			if (getCurrentLayer_().features.size()==0)
			{
				layers_.resize(getLayerCount()-1);
				if (current_layer_!=0) current_layer_ = current_layer_-1;
				QMessageBox::critical(this,"Error","Cannot add an empty dataset. Aborting!");
				return false;
			}
		}
		
		//overall values update
		recalculateRanges_(0,1,2);
		//cout << "New data range: " << overall_data_range_ << endl;
		
		resetZoom(false); //no repaint as this is done in intensityModeChange_() anyway
		
		if (getLayerCount()==2)
		{
			setIntensityMode(IM_PERCENTAGE);
		}
		intensityModeChange_();

		emit layerActivated(this);

		//set watch on the file
		if (File::exists(getCurrentLayer().filename))
		{
			watcher_->addFile(getCurrentLayer().filename.toQString());
		}

		return true;
	}
	
	void Spectrum2DCanvas::removeLayer(int layer_index )
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()))
		{
			return;
		}

		//unselect all peaks
		selected_peak_.clear();
		measurement_start_.clear();
	
		//remove the data
		layers_.erase(layers_.begin()+layer_index);
		
		//update visible area and boundaries
		recalculateRanges_(0,1,2);
		
		resetZoom(false); //no repaint as this is done in intensityModeChange_() anyway

		//update current layer if it became invalid
		if (current_layer_!=0 && current_layer_ >= getLayerCount()) current_layer_ = getLayerCount()-1;

		if (layers_.empty())
		{
			overall_data_range_ = DRange<3>::empty;
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
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

		//unselect all peaks
		selected_peak_.clear();
		measurement_start_.clear();

		current_layer_ = layer_index;
		emit layerActivated(this);
		
		update_(__PRETTY_FUNCTION__);
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
						for (ExperimentType::ConstAreaIterator it = getLayer(i).peaks.areaBeginConst(visible_area_.min()[1],visible_area_.max()[1],visible_area_.min()[0],visible_area_.max()[0]); 
								 it != getLayer(i).peaks.areaEndConst(); 
								 ++it)
						{
							PeakIndex pi = it.getPeakIndex();
							if (it->getIntensity() > local_max && getLayer(i).filters.passes(getLayer(i).peaks[pi.spectrum],pi.peak))
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
									 getLayer(i).filters.passes(*it) &&
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
		//Only fill background if no layer is present
		if (getLayerCount()==0)
		{
			QPainter painter;
			painter.begin(this);
			painter.fillRect(0,0,this->width(),this->height(),QColor(param_.getValue("background_color").toQString()));
			painter.end();
			e->accept();
			return;
		}

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
						if (getLayerFlag(i,LayerData::F_HULLS))
						{
							paintTraceConvexHulls_(i, painter);
						}
						if (getLayerFlag(i,LayerData::F_HULL))
						{
							paintFeatureConvexHulls_(i, painter);
						}
						paintDots_(i, painter);
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
		if (action_mode_==AM_MEASURE && measurement_start_.isValid())
		{
			painter.setPen(Qt::black);
			
			QPoint line_begin, line_end;
			
			if (selected_peak_.isValid())
			{
				if (getCurrentLayer().type!=LayerData::DT_PEAK)
				{
					dataToWidget_(selected_peak_.getFeature(getCurrentLayer().features).getMZ(),selected_peak_.getFeature(getCurrentLayer().features).getRT(),line_begin);
				}
				else
				{
					dataToWidget_(selected_peak_.getPeak(getCurrentLayer().peaks).getMZ(),selected_peak_.getSpectrum(getCurrentLayer().peaks).getRT(),line_begin);
				}
			}
			else
			{
				line_begin = last_mouse_pos_;
			}
			if (getCurrentLayer().type!=LayerData::DT_PEAK)
			{
				dataToWidget_(measurement_start_.getFeature(getCurrentLayer().features).getMZ(),measurement_start_.getFeature(getCurrentLayer().features).getRT(),line_end);
			}
			else
			{
				dataToWidget_(measurement_start_.getPeak(getCurrentLayer().peaks).getMZ(),measurement_start_.getSpectrum(getCurrentLayer().peaks).getRT(),line_end);
			}
			painter.drawLine(line_begin, line_end);
			
			highlightPeak_(painter, measurement_start_);
		}
		
		//draw selected peak
		if (action_mode_==AM_MEASURE || action_mode_==AM_TRANSLATE) highlightPeak_(painter, selected_peak_);
		
		//draw convex hull of selected peak
		if (selected_peak_.isValid() && getCurrentLayer().type!=LayerData::DT_PEAK)
		{
			painter.setPen(QPen(Qt::red, 2));

			paintConvexHulls_(selected_peak_.getFeature(getCurrentLayer().features).getConvexHulls(),painter);
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

		if (e->button() == Qt::LeftButton)
		{
			if (action_mode_ == AM_MEASURE)
			{
				if (selected_peak_.isValid())
				{
					measurement_start_ = selected_peak_;
				}
				else
				{
					measurement_start_.clear();
				}
			}
			else if (action_mode_ == AM_ZOOM)
			{
				rubber_band_.setGeometry(e->pos().x(),e->pos().y(),0,0);
				rubber_band_.show();
			}
		}
	}

	void Spectrum2DCanvas::mouseMoveEvent(QMouseEvent* e)
	{
		QPoint pos = e->pos();
	  PeakIndex near_peak = findNearestPeak_(pos);

		//highlight current peak and display peak coordinates
		if (action_mode_==AM_MEASURE || (action_mode_==AM_TRANSLATE && !e->buttons() & Qt::LeftButton))
		{
			//highlight peak
			selected_peak_ = near_peak; 
			update_(__PRETTY_FUNCTION__);
			
			//show coordinates
			if (selected_peak_.isValid())
			{
				if (getCurrentLayer().type!=LayerData::DT_PEAK)
				{
					//coordinates
					const FeatureMapType::FeatureType& f = selected_peak_.getFeature(getCurrentLayer().features);
					emit sendCursorStatus(f.getMZ(), f.getIntensity(), f.getRT());
					//additional feature info
					String status;
					status = status + "Quality: " + f.getOverallQuality();
					if (f.getCharge()!=0)
					{
						status = status + " Charge: " + f.getCharge();
					}
					//add meta info
					std::vector<String> keys;
					f.getKeys(keys);
					for (UInt m=0; m<keys.size(); ++m)
					{
						status = status + " " + keys[m] + ": " + (String)(f.getMetaValue(keys[m]));
					}
					emit sendStatusMessage(status, 0);
				}
				else
				{
					//coordinates
					const ExperimentType::PeakType& p = selected_peak_.getPeak(getCurrentLayer().peaks);
					const ExperimentType::SpectrumType& s = selected_peak_.getSpectrum(getCurrentLayer().peaks);
					emit sendCursorStatus(p.getMZ(), p.getIntensity(), s.getRT());
					//meta info
					String status;
					for (UInt m=0; m<s.getMetaDataArrays().size();++m)
					{
						status += s.getMetaDataArrays()[m].getName() + ": " + s.getMetaDataArrays()[m][selected_peak_.peak] + " ";
					}
					emit sendStatusMessage(status, 0);
				}
			}
		}
		else
		{
			//Zoom mode => no peak should be select
			selected_peak_.clear(); 
			update_(__PRETTY_FUNCTION__);
		}

		//show current mouse cordinates
		if (action_mode_==AM_ZOOM || !near_peak.isValid())
		{
			PointType pnt = widgetToData_(pos); 	 
			emit sendCursorStatus( pnt[0], -1.0, pnt[1]); 	 
		}

		if (action_mode_==AM_MEASURE)
		{
			last_mouse_pos_ = pos;
			//show coordinates of nearby peak or current position
			if (near_peak.isValid()) // a peak is nearby
			{
				//if a valid range is selected, show the differences
				if ((e->buttons() & Qt::LeftButton) && measurement_start_.isValid())
				{
					if (getCurrentLayer().type!=LayerData::DT_PEAK)
					{
						const FeatureMapType::FeatureType& f1 = measurement_start_.getFeature(getCurrentLayer().features);
						const FeatureMapType::FeatureType& f2 = selected_peak_.getFeature(getCurrentLayer().features);
						emit sendStatusMessage(QString("Measured: dRT = %1, dMZ = %3, Intensity ratio = %2").arg(f2.getRT()-f1.getRT()).arg(f2.getIntensity()/f1.getIntensity()).arg(f2.getMZ()-f1.getMZ()).toStdString(), 0);
					}
					else
					{
						const ExperimentType::PeakType& p1 = measurement_start_.getPeak(getCurrentLayer().peaks);
						const ExperimentType::SpectrumType& s1 = measurement_start_.getSpectrum(getCurrentLayer().peaks);
						const ExperimentType::PeakType& p2 = selected_peak_.getPeak(getCurrentLayer().peaks);
						const ExperimentType::SpectrumType& s2 = selected_peak_.getSpectrum(getCurrentLayer().peaks);
						emit sendStatusMessage(QString("Measured: dRT = %1, dMZ = %3, Intensity ratio = %2").arg(s2.getRT()-s1.getRT()).arg(p2.getIntensity()/p1.getIntensity()).arg(p2.getMZ()-p1.getMZ()).toStdString(), 0);
					}
				}
			}
		}
    else if (action_mode_ == AM_ZOOM)
		{
			//if mouse button is held down, enlarge the selection
			if (e->buttons() & Qt::LeftButton)
			{
				rubber_band_.setGeometry(last_mouse_pos_.x(), last_mouse_pos_.y(), pos.x() - last_mouse_pos_.x(), pos.y() - last_mouse_pos_.y());
				update_(__PRETTY_FUNCTION__);
			}
		}
		else if (action_mode_ == AM_TRANSLATE)
		{
			if (e->buttons() & Qt::LeftButton)
			{
				//calculate data coordinates of shift
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
		}
	}
	
	void Spectrum2DCanvas::mouseReleaseEvent(QMouseEvent* e)
	{
		QPoint pos = e->pos();
		if (e->button() == Qt::LeftButton)
		{
			if (action_mode_==AM_MEASURE)
			{
				if (!selected_peak_.isValid())
				{
					measurement_start_.clear();
				}
				update_(__PRETTY_FUNCTION__);
				if (measurement_start_.isValid())
				{
					if (getCurrentLayer().type!=LayerData::DT_PEAK)
					{
						const FeatureMapType::FeatureType& f1 = measurement_start_.getFeature(getCurrentLayer().features);
						const FeatureMapType::FeatureType& f2 = selected_peak_.getFeature(getCurrentLayer().features);
						emit sendStatusMessage(QString("Measured: dRT = %1, dMZ = %3, Intensity ratio = %2").arg(f2.getRT()-f1.getRT()).arg(f2.getIntensity()/f1.getIntensity()).arg(f2.getMZ()-f1.getMZ()).toStdString(), 0);
					}
					else
					{
						const ExperimentType::PeakType& p1 = measurement_start_.getPeak(getCurrentLayer().peaks);
						const ExperimentType::SpectrumType& s1 = measurement_start_.getSpectrum(getCurrentLayer().peaks);
						const ExperimentType::PeakType& p2 = selected_peak_.getPeak(getCurrentLayer().peaks);
						const ExperimentType::SpectrumType& s2 = selected_peak_.getSpectrum(getCurrentLayer().peaks);
						emit sendStatusMessage(QString("Measured: dRT = %1, dMZ = %3, Intensity ratio = %2").arg(s2.getRT()-s1.getRT()).arg(p2.getIntensity()/p1.getIntensity()).arg(p2.getMZ()-p1.getMZ()).toStdString(), 0);
					}
				}
				measurement_start_.clear();
			}
		  else if (action_mode_ == AM_ZOOM)
			{
				rubber_band_.hide();
				QRect rect = rubber_band_.geometry();
				if (rect.width()!=0 && rect.height()!=0)
				{
					AreaType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
					//cout << __PRETTY_FUNCTION__ << endl;
					changeVisibleArea_(area, true, true);
				}
			}
		}
	}
	
	void Spectrum2DCanvas::contextMenuEvent(QContextMenuEvent* e)
	{
		//Abort of there are no layers
		if (layers_.empty()) return;
		
		DoubleReal rt = widgetToData_(e->pos())[1];
		
		const LayerData& layer = getCurrentLayer();

		QMenu* context_menu = new QMenu(this);
		QAction* a = 0;
		QAction* result = 0;

		//Display name and warn if current layer invisible
		String layer_name = String("Layer: ") + layer.name;
		if (!layer.visible)
		{
			layer_name += " (invisible)";
		}
		context_menu->addAction(layer_name.toQString())->setEnabled(false);
		context_menu->addSeparator();
		
		context_menu->addAction("Layer meta data");
		
		QMenu* settings_menu = new QMenu("Settings");
 		settings_menu->addAction("Show/hide grid lines");
 		settings_menu->addAction("Show/hide axis legends");
 		settings_menu->addAction("Show/hide projections")->setEnabled(layer.type==LayerData::DT_PEAK);;
 		settings_menu->addAction("Show/hide MS/MS precursors")->setEnabled(layer.type==LayerData::DT_PEAK);
		settings_menu->addSeparator();
 		settings_menu->addAction("Preferences");
		
		QMenu* save_menu = new QMenu("Save");
		save_menu->addAction("Layer");
		save_menu->addAction("Visible layer data");
		
		//-------------------PEAKS----------------------------------
		if (layer.type==LayerData::DT_PEAK)
		{
			//add surrounding survey scans
			//find nearest survey scan
			Int size = getCurrentLayer().peaks.size();
			Int current = getCurrentLayer().peaks.RTBegin(rt)-getCurrentLayer().peaks.begin();
			Int i=0;
			while (current+i<size || current-i>=0)
			{
				if(current+i<size && getCurrentLayer().peaks[current+i].getMSLevel()==1)
				{
					current = current+i;
					break;
				}
				if(current-i>=0 && getCurrentLayer().peaks[current-i].getMSLevel()==1)
				{
					current = current-i;
					break;
				}
				++i;
			}
			//search for four scans in both directions
			vector<Int> indices;
			indices.push_back(current);
			i=1;
			while (current-i>=0 && indices.size()<5)
			{
				if (getCurrentLayer().peaks[current-i].getMSLevel()==1)
				{
					indices.push_back(current-i);
				}
				++i;
			}
			i=1;
			while (current+i<size && indices.size()<9)
			{
				if ( getCurrentLayer().peaks[current+i].getMSLevel()==1)
				{
					indices.push_back(current+i);
				}
				++i;
			}
			sort(indices.rbegin(),indices.rend());
			QMenu* ms1_scans = context_menu->addMenu("Survey scan in 1D");
			QMenu* ms1_meta = context_menu->addMenu("Survey scan meta data");
			context_menu->addSeparator();
			for(i=0; i<(Int)indices.size(); ++i)
			{
				if (indices[i]==current) ms1_scans->addSeparator();
				a = ms1_scans->addAction(QString("RT: ") + QString::number(getCurrentLayer().peaks[indices[i]].getRT()));
				a->setData(indices[i]);
				if (indices[i]==current) ms1_scans->addSeparator();

				if (indices[i]==current) ms1_meta->addSeparator();
				a = ms1_meta->addAction(QString("RT: ") + QString::number(getCurrentLayer().peaks[indices[i]].getRT()));
				a->setData(indices[i]);
				if (indices[i]==current) ms1_meta->addSeparator();
			}

			//add surrounding fragment scans
			QMenu* msn_scans = new QMenu("fragment scan in 1D");
			QMenu* msn_meta = new QMenu("fragment scan meta data");
			DPosition<2> p1 = widgetToData_(e->pos()+ QPoint(10,10));
			DPosition<2> p2 = widgetToData_(e->pos()- QPoint(10,10));
			DoubleReal rt_min = min(p1[1],p2[1]);
			DoubleReal rt_max = max(p1[1],p2[1]);
			DoubleReal mz_min = min(p1[0],p2[0]);
			DoubleReal mz_max = max(p1[0],p2[0]);
			bool item_added = false;
			for (ExperimentType::ConstIterator it=getCurrentLayer().peaks.RTBegin(rt_max); it!=getCurrentLayer().peaks.RTEnd(rt_min); --it)
			{
				DoubleReal mz = it->getPrecursorPeak().getPosition()[0];
				if (it->getMSLevel()>1 && mz>=mz_min && mz<=mz_max)
				{
					a = msn_scans->addAction(QString("RT: ") + QString::number(it->getRT()) + " mz: " + QString::number(mz));
					a->setData((int)(it-getCurrentLayer().peaks.begin()));
					a = msn_meta->addAction(QString("RT: ") + QString::number(it->getRT()) + " mz: " + QString::number(mz));
					a->setData((int)(it-getCurrentLayer().peaks.begin()));
					item_added=true;
				}
			}
			if (item_added)
			{
				context_menu->addMenu(msn_scans);
				context_menu->addMenu(msn_meta);
				context_menu->addSeparator();
			}
						
			//add settings menu
			context_menu->addSeparator(); 			
			context_menu->addMenu(save_menu);
 			context_menu->addMenu(settings_menu);

			//add external context menu
			if (context_add_)
			{
				context_menu->addSeparator();
				context_menu->addMenu(context_add_);
			}
			
			//evaluate menu
			if ((result = context_menu->exec(mapToGlobal(e->pos()))))
			{
				if (result->parent()==ms1_scans  || result->parent()==msn_scans)
				{
					emit showSpectrumAs1D(result->data().toInt());
				}
				else if (result->parent()==ms1_meta || result->parent()==msn_meta)
				{
					SpectrumType& spec = currentPeakData_()[result->data().toInt()];
					MSMetaDataExplorer dlg(true, this);
		      dlg.setWindowTitle("View/Edit meta data");
		    	dlg.visualize(static_cast<SpectrumSettings&>(spec));
		      //Add MetaInfoDescriptions
		      for (UInt i=0; i<spec.getMetaDataArrays().size();++i)
		      {
		      	dlg.visualize(static_cast<MetaInfoDescription&>(spec.getMetaDataArrays()[i]));
		      }
		      dlg.exec();
				}
			}	
		}
		//-------------------FEATURES----------------------------------
		else if (layer.type==LayerData::DT_FEATURE)
		{
			//search for nearby features
			DPosition<2> p1 = widgetToData_(e->pos()+ QPoint(10,10));
			DPosition<2> p2 = widgetToData_(e->pos()- QPoint(10,10));
			DoubleReal rt_min = min(p1[1],p2[1]);
			DoubleReal rt_max = max(p1[1],p2[1]);
			DoubleReal mz_min = min(p1[0],p2[0]);
			DoubleReal mz_max = max(p1[0],p2[0]);
			
			QMenu* meta = new QMenu("Feature meta data");
			bool present = false;
			FeatureMapType& features = getCurrentLayer_().features;
			for (FeatureMapType::Iterator it = features.begin(); it!=features.end(); ++it)
			{
				if (it->getMZ() <= mz_max && it->getMZ() >= mz_min && it->getRT() <= rt_max && it->getRT() >= rt_min)
				{
					present=true;
					a = meta->addAction(QString("RT: ") + QString::number(it->getRT()) + "  m/z:" + QString::number(it->getMZ()) + "  charge:" + QString::number(it->getCharge()));
					a->setData((int)(it-features.begin()));
				}  
			}
			if (present)
			{
				context_menu->addMenu(meta);
				context_menu->addSeparator(); 			
			}
			
			//add settings menu
			context_menu->addMenu(save_menu);
 			context_menu->addMenu(settings_menu);
			

			//add external context menu
			if (context_add_)
			{
				context_menu->addSeparator();
				context_menu->addMenu(context_add_);
			}

			//evaluate menu			
			if ((result = context_menu->exec(mapToGlobal(e->pos()))))
			{
				if (result->text().left(3)=="RT:")
				{
					MSMetaDataExplorer dlg(true, this);
		      dlg.setWindowTitle("View/Edit meta data");
		      dlg.visualize(features[result->data().toInt()]);
		      
		      vector<PeptideIdentification>& ids = features[result->data().toInt()].getPeptideIdentifications();
					for (vector<PeptideIdentification>::iterator it=ids.begin(); it!=ids.end(); ++it)
					{
		    		dlg.visualize(*it);
					}
					
		      dlg.exec();
				}
			}
		}
		
		//common actions of peaks and features
		if (result)
		{
			if (result->text() == "Preferences")
			{
				showCurrentLayerPreferences();
			}
			else if (result->text() == "Show/hide grid lines")
			{
				showGridLines(!gridLinesShown());
			} 
			else if (result->text() == "Show/hide axis legends")
			{
				emit changeLegendVisibility();
			}
			else if (result->text()=="Layer" || result->text()=="Visible layer data")
			{
				saveCurrentLayer(result->text()=="Visible layer data");
			}
			else if (result->text()=="Show/hide projections")
			{
				emit toggleProjections();
			}
			else if (result->text()=="Show/hide MS/MS precursors")
			{
				setLayerFlag(LayerData::P_PRECURSORS,!getLayerFlag(LayerData::P_PRECURSORS));
			}
			else if (result->text()=="Layer meta data")
			{
				showMetaData(true);
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
		QComboBox* on_file_change = dlg.findChild<QComboBox*>("on_file_change");

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
		on_file_change->setCurrentIndex(on_file_change->findText(param_.getValue("on_file_change").toQString()));	
		
		if (dlg.exec())
		{
			param_.setValue("background_color",bg_color->getColor().name());
			param_.setValue("on_file_change", on_file_change->currentText());
			if ((mapping->currentIndex()==0 && !isMzToXAxis()) || (mapping->currentIndex()==1 && isMzToXAxis()))
			{
				mzToXAxis(!isMzToXAxis());
			}
			getCurrentLayer_().param.setValue("dot:gradient",gradient->gradient().toString());
			
			currentLayerParamtersChanged_();
		}
	}

	void Spectrum2DCanvas::currentLayerParamtersChanged_()
	{
		recalculateDotGradient_(activeLayerIndex());
		
		update_buffer_ = true;	
		update_(__PRETTY_FUNCTION__);
	}

	void Spectrum2DCanvas::saveCurrentLayer(bool visible)
	{
		const LayerData& layer = getCurrentLayer();
    	
		if (layer.type==LayerData::DT_PEAK) //peak data
		{
    	QString file_name = QFileDialog::getSaveFileName(this, "Save file", param_.getValue("default_path").toQString(),"mzData files (*.mzData);;All files (*)");
			if (!file_name.isEmpty())
			{
				//set up file adapter
				MzDataFile f;
				f.setLogType(ProgressLogger::GUI);
			
	    	if (visible) //only visible data
	    	{
					ExperimentType out;
					getVisiblePeakData(out);
					f.store(file_name,out);
				}
				else //all data
				{
					f.store(file_name,getCurrentLayer().peaks);
				}
			}
	  }
	  else //features
	  {
			QString file_name = QFileDialog::getSaveFileName(this, "Save file", param_.getValue("default_path").toQString(),"featureXML files (*.featureXML);;All files (*)");
			if (!file_name.isEmpty())
			{
		  	if (visible) //only visible data
		  	{
					FeatureMapType out;
					getVisibleFeatureData(out);
					FeatureXMLFile().store(file_name,out);
				}
				else //all data
				{
					FeatureXMLFile().store(file_name,getCurrentLayer().features);
				}
			}
	  }
	}

	void Spectrum2DCanvas::updateLayer_(UInt i)
	{
		LayerData& layer = getLayer_(i);
		
		if (layers_.back().type==LayerData::DT_PEAK) //peak data
		{
			try
			{
				FileHandler().loadExperiment(layer.filename,layer.peaks);
			}
			catch(Exception::BaseException& e)
			{
				QMessageBox::critical(this,"Error",(String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
				layer.peaks.clear();
			}		
			layer.peaks.sortSpectra(true);
			layer.peaks.updateRanges(1);
		}
		else //feature data
		{
			try
			{
				FileHandler().loadFeatures(layer.filename,layer.features);
			}
			catch(Exception::BaseException& e)
			{
				QMessageBox::critical(this,"Error",(String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
				layer.features.clear();
			}
			layer.features.updateRanges();
		}
		recalculateRanges_(0,1,2);
		resetZoom(false); //no repaint as this is done in intensityModeChange_() anyway
		intensityModeChange_();
	}


	void Spectrum2DCanvas::translateLeft_()
	{
		DoubleReal shift = 0.05 * visible_area_.width();
		DoubleReal newLo = visible_area_.minX() - shift;
		DoubleReal newHi = visible_area_.maxX() - shift;
		// check if we are falling out of bounds
		if (newLo < overall_data_range_.minX())
		{
			newLo = overall_data_range_.minX();
			newHi = newLo + visible_area_.width();
		}
		//change visible area
		changeVisibleArea_(AreaType(newLo,visible_area_.minY(),newHi,visible_area_.maxY()));
	}
	
	void Spectrum2DCanvas::translateRight_()
	{
		DoubleReal shift = 0.05 * visible_area_.width();
		DoubleReal newLo = visible_area_.minX() + shift;
		DoubleReal newHi = visible_area_.maxX() + shift;
		// check if we are falling out of bounds
		if (newHi > overall_data_range_.maxX())
		{
			newHi = overall_data_range_.maxX();
			newLo = newHi - visible_area_.width();
		}
		//change visible area
		changeVisibleArea_(AreaType(newLo,visible_area_.minY(),newHi,visible_area_.maxY()));
	}
	
	void Spectrum2DCanvas::translateForward_()
	{
		DoubleReal shift = 0.05 * visible_area_.height();
		DoubleReal newLo = visible_area_.minY() + shift;
		DoubleReal newHi = visible_area_.maxY() + shift;
		// check if we are falling out of bounds
		if (newHi > overall_data_range_.maxY())
		{
			newHi = overall_data_range_.maxY();
			newLo = newHi - visible_area_.height();
		}
		//change visible area
		changeVisibleArea_(AreaType(visible_area_.minX(),newLo,visible_area_.maxX(),newHi));
	}
	
	void Spectrum2DCanvas::translateBackward_()
	{
		DoubleReal shift = 0.05 * visible_area_.height();
		DoubleReal newLo = visible_area_.minY() - shift;
		DoubleReal newHi = visible_area_.maxY() - shift;
		// check if we are falling out of bounds
		if (newLo < overall_data_range_.minY())
		{
			newLo = overall_data_range_.minY();
			newHi = newLo + visible_area_.height();
		}
		//change visible area
		changeVisibleArea_(AreaType(visible_area_.minX(),newLo,visible_area_.maxX(),newHi));		
	}


} //namespace OpenMS

