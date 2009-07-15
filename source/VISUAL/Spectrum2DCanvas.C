// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DPrefDialog.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/DIALOGS/FeatureEditDialog.h>
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
		: SpectrumCanvas(preferences, parent)
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
		
		//connect preferences change to the right slot
		connect(this,SIGNAL(preferencesChange()),this,SLOT(currentLayerParamtersChanged_()));
	}

	Spectrum2DCanvas::~Spectrum2DCanvas()
	{
	}

	void Spectrum2DCanvas::highlightPeak_(QPainter& painter, const PeakIndex& peak)
	{
		if (!peak.isValid()) return;
		
		//determine coordinates;
		QPoint pos;
		if (getCurrentLayer().type==LayerData::DT_FEATURE)
		{
			dataToWidget_(peak.getFeature(getCurrentLayer().features).getMZ(),  peak.getFeature(getCurrentLayer().features).getRT(), pos);
		}
		else if (getCurrentLayer().type==LayerData::DT_PEAK)
		{
			dataToWidget_(peak.getPeak(getCurrentLayer().peaks).getMZ(), peak.getSpectrum(getCurrentLayer().peaks).getRT(), pos);
		}
		else
		{
			dataToWidget_(peak.getFeature(getCurrentLayer().consensus).getMZ(), peak.getFeature(getCurrentLayer().consensus).getRT(), pos);
		}
		
		//paint highlighed peak
		painter.save();
		painter.setPen(QPen(Qt::red, 2));
		painter.drawEllipse(pos.x() - 5, pos.y() - 5, 10, 10);
		
		//restore painter
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
	 	else if (getCurrentLayer().type==LayerData::DT_FEATURE)
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
		else
		{
			for (ConsensusMapType::ConstIterator i = getCurrentLayer().consensus.begin();
				   i != getCurrentLayer().consensus.end();
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

						return PeakIndex(i-getCurrentLayer().consensus.begin());
					}
				}
			}
		}
		return PeakIndex();
	}

	void Spectrum2DCanvas::paintDots_(Size layer_index, QPainter& painter)
	{
		const LayerData& layer = getLayer(layer_index);
		
		//update factors (snap and percentage)
		DoubleReal snap_factor = snap_factors_[layer_index];
		percentage_factor_ = 1.0;
		if (intensity_mode_ == IM_PERCENTAGE)
		{
			if (layer.type == LayerData::DT_PEAK && layer.peaks.getMaxInt()>0.0)
			{
				percentage_factor_ = overall_data_range_.max()[2]/layer.peaks.getMaxInt();
			}
			else if (layer.type != LayerData::DT_FEATURE && layer.features.getMaxInt()>0.0)
			{
				percentage_factor_ = overall_data_range_.max()[2]/layer.features.getMaxInt();
			}
			else if (layer.type != LayerData::DT_CONSENSUS && layer.consensus.getMaxInt()>0.0)
			{
				percentage_factor_ = overall_data_range_.max()[2]/layer.consensus.getMaxInt();
			}
		}

		//set painter to black (we operate directly on the pixels for all colored data)
		painter.setPen(Qt::black);

		//temporary variables
		Int image_width = buffer_.width();
		Int image_height = buffer_.height();
		
		if (layer.type==LayerData::DT_PEAK) //peaks
		{
			//renaming some values for readability
			const ExperimentType& map = layer.peaks;
			DoubleReal rt_min = visible_area_.min()[1];
			DoubleReal rt_max = visible_area_.max()[1];
			DoubleReal mz_min = visible_area_.min()[0];
			DoubleReal mz_max = visible_area_.max()[0];

			//determine number of pixels for each dimension
			Int rt_pixel_count = image_height;
			Int mz_pixel_count = image_width;
			if(!isMzToXAxis())
			{
				rt_pixel_count = image_width;
				mz_pixel_count = image_height;
			}
			
			//-----------------------------------------------------------------------------------------------
			//determine if we want to draw dots or crosses (this also influences the way the data is painted)
			//determine number of shown scans
			UInt scans = 0;
			for (ExperimentType::ConstIterator it = map.RTBegin(rt_min); it!=map.RTEnd(rt_max); ++it)
			{
				if (it->getMSLevel()==1) ++scans;
			}
			//estimate the number of shown peaks from the central MS1 scan of the RT range
			Int peaks = 0;
			{
				ExperimentType::ConstIterator it = map.RTBegin(rt_min) + scans/2;
				while (it!=map.end() && it->getMSLevel()!=1)
				{
					++it;
				}
				if (it!=map.end())
				{
					for  (ExperimentType::SpectrumType::ConstIterator it2 = it->MZBegin(mz_min); it2!=it->MZEnd(mz_max); ++it2)
					{
						if (layer.filters.passes(*it,it2-it->begin())) ++peaks;
					}
				}
			}
			
			//-----------------------------------------------------------------------------------------------
			//paint dots (many data points): we paint the maximum shown intensity per pixel
			if (peaks*scans>0.25*mz_pixel_count*rt_pixel_count)
			{
				//calculate pixel size in data coordinates
				DoubleReal rt_step_size = (rt_max - rt_min) / rt_pixel_count;
				DoubleReal mz_step_size = (mz_max - mz_min) / mz_pixel_count;
				
				//iterate over all pixels (RT dimension)
				Size scan_index = 0;
				for (Int rt=0; rt<rt_pixel_count; ++rt)
				{
					DoubleReal rt_start = rt_min + rt_step_size * rt; 
					DoubleReal rt_end = rt_start + rt_step_size; 
					//cout << "rt: " << rt << " (" << rt_start << " - " << rt_end << ")" << endl;
					
					//determine the relevant spectra and reserve an array for the peak indices
					vector<Size> scan_indices, peak_indices;
					for (Size i=scan_index; i<map.size(); ++i)
					{
						if (map[i].getRT()>=rt_end)
						{
							scan_index = i; //store last scan index for next RT pixel
							break;
						}
						if (map[i].getMSLevel()==1 && map[i].size()>0)
						{
							scan_indices.push_back(i);
							peak_indices.push_back(map[i].MZBegin(mz_min) - map[i].begin());
						}
						//set the scan index past the end. Otherwise the last scan will be repeated for all following RTs
						if (i==map.size()-1) scan_index=i+1;
					}
					//cout << "  scans: " << scan_indices.size() << endl;

					if (scan_indices.size()==0) continue;
					
					//iterate over all pixels (m/z dimension)
					for (Int mz=0; mz<mz_pixel_count; ++mz)
					{
						DoubleReal mz_start = mz_min + mz_step_size * mz;
						DoubleReal mz_end = mz_start + mz_step_size; 
						
						//iterate over all relevant peaks in all relevant scans
						Real max = -1.0;
						for (Size i=0; i<scan_indices.size(); ++i)
						{
							Size s = scan_indices[i];
							Size p = peak_indices[i];
							for(; p<map[s].size(); ++p)
							{
								if (map[s][p].getMZ()>=mz_end) break;
								if (map[s][p].getIntensity() > max && layer.filters.passes(map[s],p))
								{
									max = map[s][p].getIntensity();
								}
							}
							peak_indices[i] = p; //store last peak index for next m/z pixel
						}
						
						//draw to buffer
						if (max>=0.0)
						{
							QPoint pos;
							dataToWidget_(mz_start + 0.5 * mz_step_size, rt_start + 0.5 * rt_step_size, pos);
							if (pos.y()<image_height && pos.x()<image_width)
							{
								buffer_.setPixel(pos.x() , pos.y(), heightColor_(max, layer.gradient, snap_factor));
							}
						}
					}
				}
			}
			//-----------------------------------------------------------------------------------------------
			//paint crosses (few data points): we paint the data on the screen
			else
			{
				for (ExperimentType::ConstAreaIterator i = map.areaBeginConst(rt_min,rt_max,mz_min,mz_max);
						 i != map.areaEndConst();
						 ++i)
				{
					PeakIndex pi = i.getPeakIndex();
					if (layer.filters.passes(map[pi.spectrum],pi.peak))
					{
						QPoint pos;
						dataToWidget_(i->getMZ(), i.getRT(), pos);
						if (pos.x()>0 && pos.y()>0 && pos.x()<image_width-1 && pos.y()<image_height-1)
						{
							QRgb color = heightColor_(i->getIntensity(), layer.gradient, snap_factor);
							buffer_.setPixel(pos.x() ,pos.y() ,color);
							buffer_.setPixel(pos.x()-1 ,pos.y() ,color);
							buffer_.setPixel(pos.x()+1 ,pos.y() ,color);
							buffer_.setPixel(pos.x() ,pos.y()-1 ,color);
							buffer_.setPixel(pos.x() ,pos.y()+1 ,color);
						}
					}
				}
			}

			//-----------------------------------------------------------------
			//draw precursor peaks
			if (getLayerFlag(layer_index,LayerData::P_PRECURSORS))
			{
				for (ExperimentType::ConstIterator i = map.RTBegin(visible_area_.min()[1]);
						 i != map.RTEnd(visible_area_.max()[1]);
						 ++i)
				{
					//this is an MS/MS scan
					if (i->getMSLevel()==2 && !i->getPrecursors().empty())
					{
						ExperimentType::ConstIterator prec=map.getPrecursorSpectrum(i);
						if (prec!=map.end())
						{
							QPoint pos;
							dataToWidget_(i->getPrecursors()[0].getMZ(), prec->getRT(),pos);
							painter.drawLine(pos.x(),pos.y()+2,pos.x()+2,pos.y());
							painter.drawLine(pos.x()+2,pos.y(),pos.x(),pos.y()-2);
							painter.drawLine(pos.x(),pos.y()-2,pos.x()-2,pos.y());
							painter.drawLine(pos.x()-2,pos.y(),pos.x(),pos.y()+2);
						}
					}
				}
			}
		}
		else if (layer.type==LayerData::DT_FEATURE) //features
		{
			bool show_label = (layer.label!=LayerData::L_NONE);
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
					//determine color
					QRgb color;
					if (i->metaValueExists(5))
					{
						color = QColor(i->getMetaValue(5).toQString()).rgb();
					}
					else
					{
						color = heightColor_(i->getIntensity(), layer.gradient, snap_factor);
					}
					//paint
					QPoint pos;
					dataToWidget_(i->getMZ(),i->getRT(),pos);
					if (pos.x()>0 && pos.y()>0 && pos.x()<image_width-1 && pos.y()<image_height-1)
					{
						buffer_.setPixel(pos.x()   ,pos.y()   ,color);
						buffer_.setPixel(pos.x()-1 ,pos.y()   ,color);
						buffer_.setPixel(pos.x()+1 ,pos.y()   ,color);
						buffer_.setPixel(pos.x()   ,pos.y()-1 ,color);
						buffer_.setPixel(pos.x()   ,pos.y()+1 ,color);
					}
					//labels
					if (show_label)
					{
						if (layer.label==LayerData::L_INDEX)
						{
							painter.drawText(pos.x()+10,pos.y()+10,QString::number(num));
						}
						else if (layer.label==LayerData::L_ID)
						{
							if (i->getPeptideIdentifications().size() && i->getPeptideIdentifications()[0].getHits().size())
							{
								painter.drawText(pos.x()+10,pos.y()+10,i->getPeptideIdentifications()[0].getHits()[0].getSequence().toString().toQString());
							}
						}
						else if (layer.label==LayerData::L_META_LABEL)
						{
							painter.drawText(pos.x()+10,pos.y()+10,i->getMetaValue(3).toQString());
						}
					}
				}
				++num;
			}
		}
		else // consensus features
		{
			for (ConsensusMapType::ConstIterator i = layer.consensus.begin();
				   i != layer.consensus.end();
				   ++i)
			{
				if ( i->getRT() >= visible_area_.min()[1] &&
						 i->getRT() <= visible_area_.max()[1] &&
						 i->getMZ() >= visible_area_.min()[0] &&
						 i->getMZ() <= visible_area_.max()[0] &&
						 layer.filters.passes(*i))
				{
					//determine color
					QRgb color = heightColor_(i->getIntensity(), layer.gradient, snap_factor);
					//paint
					QPoint pos;
					dataToWidget_(i->getMZ(),i->getRT(),pos);
					if (pos.x()>0 && pos.y()>0 && pos.x()<image_width-1 && pos.y()<image_height-1)
					{
						buffer_.setPixel(pos.x()   ,pos.y()   ,color);
						buffer_.setPixel(pos.x()-1 ,pos.y()   ,color);
						buffer_.setPixel(pos.x()+1 ,pos.y()   ,color);
						buffer_.setPixel(pos.x()   ,pos.y()-1 ,color);
						buffer_.setPixel(pos.x()   ,pos.y()+1 ,color);
					}
				}
			}
		}
	}

	void Spectrum2DCanvas::paintTraceConvexHulls_(Size layer_index, QPainter& painter)
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

	void Spectrum2DCanvas::paintFeatureConvexHulls_(Size layer_index, QPainter& painter)
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
				points.resize((int)hull.getPoints().size());

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
		for (Size hull=0; hull<hulls.size(); ++hull)
		{
			points.resize((int)hulls[hull].getPoints().size());
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

	void Spectrum2DCanvas::paintConsensusElements_(Size layer_index, QPainter& p)
	{
		const LayerData& layer = getLayer(layer_index);

		for (ConsensusMapType::ConstIterator i = layer.consensus.begin(); i != layer.consensus.end(); ++i)
		{
			paintConsensusElement_(layer_index, *i, p, true);
		}
	}

	void Spectrum2DCanvas::paintConsensusElement_(Size layer_index, const ConsensusFeature& cf, QPainter& p, bool use_buffer)
	{
		Int image_width = buffer_.width();
		Int image_height = buffer_.height();

		const LayerData& layer = getLayer(layer_index);

		if ( isConsensusFeatureVisible_(cf, layer_index) && layer.filters.passes(cf))
		{
			//calculate position of consensus feature (centroid)
			QPoint consensus_pos;
			dataToWidget_(cf.getMZ(), cf.getRT(),consensus_pos);
			//iterate over elements
			for (ConsensusFeature::HandleSetType::const_iterator element=cf.begin(); element!=cf.end(); ++element)
			{
				//calculate position of consensus element
				QPoint pos;
				dataToWidget_(element->getMZ(), element->getRT(),pos);
				//paint line
				p.drawLine(consensus_pos,pos);
				//paint point
				if (pos.x()>0 && pos.y()>0 && pos.x()<image_width-1 && pos.y()<image_height-1)
				{
					// use buffer only when not highlighting
					if (use_buffer)
					{
						buffer_.setPixel(pos.x()   ,pos.y()   ,Qt::black);
						buffer_.setPixel(pos.x()-1 ,pos.y()   ,Qt::black);
						buffer_.setPixel(pos.x()+1 ,pos.y()   ,Qt::black);
						buffer_.setPixel(pos.x()   ,pos.y()-1 ,Qt::black);
						buffer_.setPixel(pos.x()   ,pos.y()+1 ,Qt::black);
					}
					else
					{
						p.drawPoint(pos.x()   ,pos.y());
						p.drawPoint(pos.x()-1 ,pos.y());
						p.drawPoint(pos.x()+1 ,pos.y());
						p.drawPoint(pos.x()   ,pos.y()-1);
						p.drawPoint(pos.x()   ,pos.y()+1);
					}
				}
			}
		}

	}


	bool Spectrum2DCanvas::isConsensusFeatureVisible_(const ConsensusFeature& ce, Size layer_index)
	{
		// check the centroid first
		if (ce.getRT() >= visible_area_.min()[1] &&
				ce.getRT() <= visible_area_.max()[1] &&
				ce.getMZ() >= visible_area_.min()[0] &&
				ce.getMZ() <= visible_area_.max()[0])
		{
			return true;
		}

		// if element-flag is set, check if any of the consensus elements is visible
		if (getLayerFlag(layer_index, LayerData::C_ELEMENTS))
		{
			ConsensusFeature::HandleSetType::const_iterator element=ce.getFeatures().begin();
			for (; element != ce.getFeatures().end(); ++element)
			{
				if (element->getRT() >= visible_area_.min()[1] &&
						element->getRT() <= visible_area_.max()[1] &&
						element->getMZ() >= visible_area_.min()[0] &&
						element->getMZ() <= visible_area_.max()[0])
				{
					return true;
				}
			}
		}
		return false;
	}

	void Spectrum2DCanvas::intensityModeChange_()
	{
		for (Size i=0; i<layers_.size();++i)
		{
			recalculateDotGradient_(i);
		}
		SpectrumCanvas::intensityModeChange_();
	}

	void Spectrum2DCanvas::recalculateDotGradient_(Size layer)
	{
		getLayer_(layer).gradient.fromString(getLayer_(layer).param.getValue("dot:gradient"));
		getLayer_(layer).gradient.activatePrecalculationMode(getMinIntensity(layer), overall_data_range_.max()[2], param_.getValue("interpolation_steps"));
	}

	void Spectrum2DCanvas::updateProjections()
	{
		//try to find the right (peak) layer to project
		const LayerData* layer = &(getCurrentLayer()); //first, try current layer
		if (layer->type != LayerData::DT_PEAK) //second, check if more than one peak layer is present
		{
			UInt peak_layer_count = 0;
			Size last_peak_layer = 0;
			for (Size i=0; i<getLayerCount(); ++i)
			{
				if (getLayer(i).type==LayerData::DT_PEAK)
				{
					++peak_layer_count;
					last_peak_layer = i;
				}
			}
			if (peak_layer_count==1)
			{
				layer = &(getLayer(last_peak_layer));
			}
		}
		if (layer->type != LayerData::DT_PEAK) //third, check if more than one peak layer is visible
		{
			UInt peak_layer_count = 0;
			Size last_peak_layer = 0;
			for (Size i=0; i<getLayerCount(); ++i)
			{
				if (getLayer(i).type==LayerData::DT_PEAK && getLayer(i).visible)
				{
					++peak_layer_count;
					last_peak_layer = i;
				}
			}
			if (peak_layer_count==1)
			{
				layer = &(getLayer(last_peak_layer));
			}
		}
		if (layer->type != LayerData::DT_PEAK)//forth, abort of no peak layer or several peak layers are present
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

	  float prec = 0.05f;
		float mult = 1.0/prec;


		for (ExperimentType::ConstAreaIterator i = layer->peaks.areaBeginConst(visible_area_.min()[1],visible_area_.max()[1],visible_area_.min()[0],visible_area_.max()[0]);
				 i != layer->peaks.areaEndConst();
				 ++i)
		{
			PeakIndex pi = i.getPeakIndex();
			if (layer->filters.passes(layer->peaks[pi.spectrum],pi.peak))
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

		Size i = 2;
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
		else if (layers_.back().type==LayerData::DT_FEATURE)//feature data
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
		else //consensus feature data
		{
			getCurrentLayer_().consensus.updateRanges();

			//Abort if no data points are contained
			if (getCurrentLayer_().consensus.size()==0)
			{
				layers_.resize(getLayerCount()-1);
				if (current_layer_!=0) current_layer_ = current_layer_-1;
				QMessageBox::critical(this,"Error","Cannot add an empty dataset. Aborting!");
				return false;
			}
		}

		//Warn if negative intensities are contained
		if (getMinIntensity(current_layer_)<0.0)
		{
			QMessageBox::warning(this,"Warning","This dataset contains negative intensities. Use it at your own risk!");
		}

		//overall values update
		recalculateRanges_(0,1,2);
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

	void Spectrum2DCanvas::removeLayer(Size layer_index )
	{
		if (layer_index >= getLayerCount())
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
	void Spectrum2DCanvas::activateLayer(Size layer_index)
	{
		if (layer_index >= getLayerCount() || layer_index==current_layer_)
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
		snap_factors_ = vector<DoubleReal>(getLayerCount(),1.0);

		if (intensity_mode_ == IM_SNAP)
		{
			for (Size i=0; i<getLayerCount(); i++)
			{
				if (getLayer(i).visible)
				{
					DoubleReal local_max  = -numeric_limits<DoubleReal>::max();
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
					else if (getLayer(i).type==LayerData::DT_FEATURE) //features
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
					else
					{
						for (ConsensusMapType::ConstIterator it = getLayer(i).consensus.begin();
							   it != getLayer(i).consensus.end();
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
					if (local_max>0.0)
					{
						snap_factors_[i] = overall_data_range_.max()[2]/local_max;
					}
				}
			}
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
	  cout << "  Visible area -- m/z: " << visible_area_.minX() << " - " << visible_area_.maxX() << " rt: " << visible_area_.minY() << " - " << visible_area_.maxY() << endl;
	  cout << "  Overall area -- m/z: " << overall_data_range_.min()[0] << " - " << overall_data_range_.max()[0] << " rt: " << overall_data_range_.min()[1] << " - " << overall_data_range_.max()[1] << endl;
#endif
		
		//timing
		QTime overall_timer;
		if (show_timing_)
		{
			overall_timer.start();
			if (update_buffer_)
			{
				cout << "Updating buffer:" << endl;				
			}
			else
			{
				cout << "Copying buffer:" << endl;				
			}
		}
		
		QPainter painter;
		if (update_buffer_)
		{
			update_buffer_ = false;

			//recalculate snap factor
			recalculateSnapFactor_();

			buffer_.fill(QColor(param_.getValue("background_color").toQString()).rgb());
			painter.begin(&buffer_);
			QTime layer_timer;
			
			for (Size i=0; i<getLayerCount(); i++)
			{
				//timing
				if (show_timing_)
				{
					layer_timer.start();
				}
				
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
					else
					{
						if (getLayerFlag(i,LayerData::C_ELEMENTS))
						{
							paintConsensusElements_(i, painter);
						}
						paintDots_(i, painter);
					}
				}
				//timing
				if (show_timing_)
				{
					cout << "  -layer " << i << " time: " << layer_timer.elapsed() << " ms" << endl;
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
			painter.drawImage(rects[i].topLeft(), buffer_, rects[i]);
		}

		//draw measurement peak
		if (action_mode_==AM_MEASURE && measurement_start_.isValid())
		{
			painter.setPen(Qt::black);

			QPoint line_begin, line_end;

			if (selected_peak_.isValid())
			{
				if (getCurrentLayer().type==LayerData::DT_FEATURE)
				{
					dataToWidget_(selected_peak_.getFeature(getCurrentLayer().features).getMZ(),selected_peak_.getFeature(getCurrentLayer().features).getRT(),line_begin);
				}
				else if (getCurrentLayer().type==LayerData::DT_PEAK)
				{
					dataToWidget_(selected_peak_.getPeak(getCurrentLayer().peaks).getMZ(),selected_peak_.getSpectrum(getCurrentLayer().peaks).getRT(),line_begin);
				}
				else
				{
					dataToWidget_(selected_peak_.getFeature(getCurrentLayer().consensus).getMZ(),selected_peak_.getFeature(getCurrentLayer().consensus).getRT(),line_begin);
				}
			}
			else
			{
				line_begin = last_mouse_pos_;
			}
			if (getCurrentLayer().type==LayerData::DT_FEATURE)
			{
				dataToWidget_(measurement_start_.getFeature(getCurrentLayer().features).getMZ(),measurement_start_.getFeature(getCurrentLayer().features).getRT(),line_end);
			}
			else if (getCurrentLayer().type==LayerData::DT_PEAK)
			{
				dataToWidget_(measurement_start_.getPeak(getCurrentLayer().peaks).getMZ(),measurement_start_.getSpectrum(getCurrentLayer().peaks).getRT(),line_end);
			}
			else
			{
				dataToWidget_(measurement_start_.getFeature(getCurrentLayer().consensus).getMZ(),measurement_start_.getFeature(getCurrentLayer().consensus).getRT(),line_end);
			}

			painter.drawLine(line_begin, line_end);

			highlightPeak_(painter, measurement_start_);
		}
		
		//draw convex hulls or consensus feature elements
		if (selected_peak_.isValid())
		{
			if (getCurrentLayer().type==LayerData::DT_FEATURE)
			{
				painter.setPen(QPen(Qt::red, 2));
				paintConvexHulls_(selected_peak_.getFeature(getCurrentLayer().features).getConvexHulls(),painter);
			}
			else if (getCurrentLayer().type==LayerData::DT_CONSENSUS && getLayerFlag(current_layer_,LayerData::C_ELEMENTS))
			{
				painter.setPen(QPen(Qt::red, 2));
				paintConsensusElement_(current_layer_, selected_peak_.getFeature(getCurrentLayer().consensus),painter,false);
			}
		}

		if (action_mode_==AM_MEASURE || action_mode_==AM_TRANSLATE)
		{
			highlightPeak_(painter, selected_peak_);
		}
		
		//draw delta for measuring
		if (action_mode_==AM_MEASURE && measurement_start_.isValid() && selected_peak_.isValid())
		{
			drawDeltas_(painter, measurement_start_, selected_peak_, true);
		}
		else
		{
			drawCoordinates_(painter, selected_peak_, true);
		}

		painter.end();
#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
		if (show_timing_)
		{
			cout << "  -overall time: " << overall_timer.elapsed() << " ms" << endl << endl;
		}
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
				//translate (if not moving features)
				if (!getCurrentLayer().type==LayerData::DT_FEATURE || !selected_peak_.isValid())
				{
					rubber_band_.setGeometry(QRect(e->pos(),QSize()));
					rubber_band_.show();
				}
			}
		}
	}

	void Spectrum2DCanvas::mouseMoveEvent(QMouseEvent* e)
	{
		grabKeyboard(); // (re-)grab keyboard after it has been released by unhandled key
		QPoint pos = e->pos();
		PointType data_pos = widgetToData_(pos);
		emit sendCursorStatus( data_pos[0], data_pos[1]);
					
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
				if (getCurrentLayer().type==LayerData::DT_FEATURE)
				{
					//coordinates
					const FeatureMapType::FeatureType& f = selected_peak_.getFeature(getCurrentLayer().features);
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
					for (Size m=0; m<keys.size(); ++m)
					{
						status = status + " " + keys[m] + ": " + (String)(f.getMetaValue(keys[m]));
					}
					emit sendStatusMessage(status, 0);
				}
				else if (getCurrentLayer().type==LayerData::DT_PEAK)
				{
					//meta info
					String status;
					const ExperimentType::SpectrumType& s = selected_peak_.getSpectrum(getCurrentLayer().peaks);
					for (Size m=0; m<s.getFloatDataArrays().size();++m)
					{
						status += s.getFloatDataArrays()[m].getName() + ": " + s.getFloatDataArrays()[m][selected_peak_.peak] + " ";
					}
					for (Size m=0; m<s.getStringDataArrays().size();++m)
					{
						status += s.getStringDataArrays()[m].getName() + ": " + s.getStringDataArrays()[m][selected_peak_.peak] + " ";
					}
					emit sendStatusMessage(status, 0);
				}
				else // ConsensusFeature
				{
					//additional feature info
					const ConsensusFeature& f = selected_peak_.getFeature(getCurrentLayer().consensus);
					String status;
					status = status + "Quality: " + f.getQuality();
					if (f.getCharge()!=0)
					{
						status = status + " Charge: " + f.getCharge();
					}
					//add meta info
					std::vector<String> keys;
					f.getKeys(keys);
					for (Size m=0; m<keys.size(); ++m)
					{
						status = status + " " + keys[m] + ": " + (String)(f.getMetaValue(keys[m]));
					}
					emit sendStatusMessage(status, 0);
				}
			}
		}
		else if (action_mode_==AM_ZOOM)
		{
			//Zoom mode => no peak should be selected
			selected_peak_.clear();
			update_(__PRETTY_FUNCTION__);
		}

		if (action_mode_==AM_MEASURE)
		{
			last_mouse_pos_ = pos;
		}
    else if (action_mode_ == AM_ZOOM)
		{
			//if mouse button is held down, enlarge the selection
			if (e->buttons() & Qt::LeftButton)
			{
				rubber_band_.setGeometry(QRect(last_mouse_pos_,pos).normalized());
				rubber_band_.show(); //if the mouse button is pressed before the zoom key is pressed

				update_(__PRETTY_FUNCTION__);
			}
		}
		else if (action_mode_ == AM_TRANSLATE)
		{
			if (e->buttons() & Qt::LeftButton)
			{
				if (getCurrentLayer().type==LayerData::DT_FEATURE && selected_peak_.isValid()) //move feature
				{
					PointType new_data = widgetToData_(pos);
					DoubleReal mz = new_data[0];
					DoubleReal rt = new_data[1];

					//restrict the movement to the data range
					mz = max(mz,overall_data_range_.min()[0]);
					mz = min(mz,overall_data_range_.max()[0]);
					rt = max(rt,overall_data_range_.min()[1]);
					rt = min(rt,overall_data_range_.max()[1]);

					getCurrentLayer_().features[selected_peak_.peak].setRT(rt);
					getCurrentLayer_().features[selected_peak_.peak].setMZ(mz);

					update_buffer_ = true;
					update_(__PRETTY_FUNCTION__);
					modificationStatus_(activeLayerIndex(), true);
				}
				else //translate
				{
					//calculate data coordinates of shift
					PointType old_data = widgetToData_(last_mouse_pos_);
					PointType new_data = widgetToData_(pos);
					//calculate x shift
					DoubleReal shift = old_data.getX() - new_data.getX();
					DoubleReal newLoX = visible_area_.minX() + shift;
					DoubleReal newHiX = visible_area_.maxX() + shift;
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
					DoubleReal newLoY = visible_area_.minY() + shift;
					DoubleReal newHiY = visible_area_.maxY() + shift;
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
				measurement_start_.clear();
				update_(__PRETTY_FUNCTION__);
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

		//-------------------PEAKS----------------------------------
		if (layer.type==LayerData::DT_PEAK)
		{
			//add settings
			settings_menu->addSeparator();
 			settings_menu->addAction("Show/hide projections");
 			settings_menu->addAction("Show/hide MS/MS precursors");

			//add surrounding survey scans
			//find nearest survey scan
			SignedSize size = getCurrentLayer().peaks.size();
			Int current = getCurrentLayer().peaks.RTBegin(rt)-getCurrentLayer().peaks.begin();
			SignedSize i=0;
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
			for (ExperimentType::ConstIterator it=getCurrentLayer().peaks.RTBegin(rt_min); it!=getCurrentLayer().peaks.RTEnd(rt_max); ++it)
			{
				DoubleReal mz = 0.0;
				if (!it->getPrecursors().empty()) mz = it->getPrecursors()[0].getMZ();
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

			finishContextMenu_(context_menu, settings_menu);

			//evaluate menu
			if ((result = context_menu->exec(mapToGlobal(e->pos()))))
			{
				if (result->parent()==ms1_scans  || result->parent()==msn_scans)
				{
					emit showSpectrumAs1D(result->data().toInt());
				}
				else if (result->parent()==ms1_meta || result->parent()==msn_meta)
				{
					showMetaData(true, result->data().toInt());
				}
			}
		}
		//-------------------FEATURES----------------------------------
		else if (layer.type==LayerData::DT_FEATURE)
		{
			//add settings
			settings_menu->addSeparator();
 			settings_menu->addAction("Show/hide convex hull");
 			settings_menu->addAction("Show/hide trace convex hulls");
 			settings_menu->addAction("Show/hide numbers/labels");

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
			//featre meta data menu
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

			finishContextMenu_(context_menu, settings_menu);

			//evaluate menu
			if ((result = context_menu->exec(mapToGlobal(e->pos()))))
			{
				if (result->text().left(3)=="RT:")
				{
					showMetaData(true,result->data().toInt());
				}
			}
		}
		//-------------------CONSENSUS FEATURES----------------------------------
		else if (layer.type==LayerData::DT_CONSENSUS)
		{
			//add settings
			settings_menu->addSeparator();
 			settings_menu->addAction("Show/hide elements");

			//search for nearby features
			DPosition<2> p1 = widgetToData_(e->pos()+ QPoint(10,10));
			DPosition<2> p2 = widgetToData_(e->pos()- QPoint(10,10));
			DoubleReal rt_min = min(p1[1],p2[1]);
			DoubleReal rt_max = max(p1[1],p2[1]);
			DoubleReal mz_min = min(p1[0],p2[0]);
			DoubleReal mz_max = max(p1[0],p2[0]);

			QMenu* consens_meta = new QMenu("Consensus meta data");
			bool present = false;
			ConsensusMapType& features = getCurrentLayer_().consensus;
			//consensus feature meta data menu
			for (ConsensusMapType::Iterator it = features.begin(); it!=features.end(); ++it)
			{
				if (it->getMZ() <= mz_max && it->getMZ() >= mz_min && it->getRT() <= rt_max && it->getRT() >= rt_min)
				{
					present=true;

					a = consens_meta->addAction(QString("RT: ") + QString::number(it->getRT()) + "  m/z:" + QString::number(it->getMZ()) + "  charge:" + QString::number(it->getCharge()));
					a->setData((int)(it-features.begin()));
				}
			}
			if (present)
			{
				context_menu->addMenu(consens_meta);
				context_menu->addSeparator();
			}

			finishContextMenu_(context_menu, settings_menu);

			if ((result = context_menu->exec(mapToGlobal(e->pos()))))
			{
				if (result->text().left(3)=="RT:")
				{
					showMetaData(true, result->data().toInt());
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
			else if (result->text()=="As image")
			{
				spectrum_widget_->saveAsImage();
			}
			else if (result->text()=="Show/hide projections")
			{
				emit toggleProjections();
			}
			else if (result->text()=="Show/hide MS/MS precursors")
			{
				setLayerFlag(LayerData::P_PRECURSORS,!getLayerFlag(LayerData::P_PRECURSORS));
			}
			else if (result->text()=="Show/hide convex hull")
			{
				setLayerFlag(LayerData::F_HULL,!getLayerFlag(LayerData::F_HULL));
			}
			else if (result->text()=="Show/hide trace convex hulls")
			{
				setLayerFlag(LayerData::F_HULLS,!getLayerFlag(LayerData::F_HULLS));
			}
			else if (result->text()=="Show/hide numbers/labels")
			{
				if (layer.label==LayerData::L_NONE) 
				{
					getCurrentLayer_().label=LayerData::L_META_LABEL;
				}
				else 
				{
					getCurrentLayer_().label=LayerData::L_NONE;
				}
			}
			else if (result->text()=="Show/hide elements")
			{
				setLayerFlag(LayerData::C_ELEMENTS,!getLayerFlag(LayerData::C_ELEMENTS));
			}
			else if (result->text()=="Layer meta data")
			{
				showMetaData(true);
			}
		}

		e->accept();
	}

	void Spectrum2DCanvas::finishContextMenu_(QMenu* context_menu, QMenu* settings_menu)
	{
		//finish settings menu
		settings_menu->addSeparator();
 		settings_menu->addAction("Preferences");

 		//create save menu
 		QMenu* save_menu = new QMenu("Save");
		save_menu->addAction("Layer");
		save_menu->addAction("Visible layer data");
		save_menu->addAction("As image");

		//add settings menu
		context_menu->addMenu(save_menu);
		context_menu->addMenu(settings_menu);

		//add external context menu
		if (context_add_)
		{
			context_menu->addSeparator();
			context_menu->addMenu(context_add_);
		}
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

		  emit preferencesChange();
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

    //determine proposed filename
		String proposed_name = param_.getValue("default_path");
    if (visible==false && layer.filename!="")
    {
    	proposed_name = layer.filename;
    }

		if (layer.type==LayerData::DT_PEAK) //peak data
		{
    	QString file_name = QFileDialog::getSaveFileName(this, "Save file", proposed_name.toQString(),"mzML files (*.mzML);;All files (*)");
			if (!file_name.isEmpty())
			{
				//set up file adapter
				MzMLFile f;
				f.setLogType(ProgressLogger::GUI);

	    	if (visible) //only visible data
	    	{
					ExperimentType out;
					getVisiblePeakData(out);
					addDataProcessing_(out, DataProcessing::FILTERING);
					f.store(file_name,out);
				}
				else //all data
				{
					f.store(file_name,layer.peaks);
				}
				modificationStatus_(activeLayerIndex(), false);
			}
	  }
	  else if (layer.type==LayerData::DT_FEATURE) //features
	  {
			QString file_name = QFileDialog::getSaveFileName(this, "Save file", proposed_name.toQString(),"FeatureXML files (*.featureXML);;All files (*)");
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
					FeatureXMLFile().store(file_name,layer.features);
				}
				modificationStatus_(activeLayerIndex(), false);
			}
	  }
	  else //consensus feature data
	  {
			QString file_name = QFileDialog::getSaveFileName(this, "Save file", proposed_name.toQString(),"ConsensusXML files (*.consensusXML);;All files (*)");
			if (!file_name.isEmpty())
			{
		  	if (visible) //only visible data
		  	{
					ConsensusMapType out;
					getVisibleConsensusData(out);
					ConsensusXMLFile().store(file_name,out);
				}
				else //all data
				{
					ConsensusXMLFile().store(file_name,layer.consensus);
				}
				modificationStatus_(activeLayerIndex(), false);
			}
	  }
	}

	void Spectrum2DCanvas::updateLayer_(Size i)
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
		else if (layers_.back().type==LayerData::DT_FEATURE) //feature data
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
		else if (layers_.back().type==LayerData::DT_CONSENSUS)  //consensus feature data
		{
			try
			{
				ConsensusXMLFile().load(layer.filename,layer.consensus);
			}
			catch(Exception::BaseException& e)
			{
				QMessageBox::critical(this,"Error",(String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
				layer.consensus.clear();
			}
			layer.consensus.updateRanges();
		}
		recalculateRanges_(0,1,2);
		resetZoom(false); //no repaint as this is done in intensityModeChange_() anyway
		intensityModeChange_();
		modificationStatus_(i, false);
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


	void Spectrum2DCanvas::keyPressEvent(QKeyEvent* e)
	{
		// CTRL+ALT+H => hidden action
		if ((e->modifiers() & Qt::ControlModifier) && (e->modifiers() & Qt::AltModifier) && (e->key()==Qt::Key_H))
		{
			/*
			//Scaling with file size (layers)
			for (UInt i=0; i<getLayerCount(); ++i)
			{
				QTime timer;
				timer.start();
				for (UInt j=0; j<10; ++j)
				{
					QPainter painter(&buffer_);
					paintDots_(i, painter);
				}
				cout << "peaks: " << getLayer(i).peaks.getSize() << " time: " << timer.elapsed() / 10.0 << endl;
			}
			
			//Scaling with resolution
			for (UInt i=250; i<3001; i+=250)
			{
				resize(i,i);
				QTime timer;
				timer.start();
				QPainter painter(&buffer_);
				paintDots_(0, painter);
				cout << "pixels: " << i << " time: " << timer.elapsed() << endl;
			}
			*/

			e->accept();
			return;
		}

		// Delete features
		LayerData& layer = getCurrentLayer_();
		if (layer.type==LayerData::DT_FEATURE && selected_peak_.isValid() && e->key()==Qt::Key_Delete)
		{
			layer.features.erase(layer.features.begin()+selected_peak_.peak);
			selected_peak_.clear();
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
			e->accept();

			modificationStatus_(activeLayerIndex(), true);
		}
		else
		{
			SpectrumCanvas::keyPressEvent(e);
		}
	}

	void Spectrum2DCanvas::keyReleaseEvent(QKeyEvent* e)
	{
		//zoom if in zoom mode and a valid rectangle is selected
		if (action_mode_==AM_ZOOM && rubber_band_.isVisible())
		{
			rubber_band_.hide();
			QRect rect = rubber_band_.geometry();
			if (rect.width()!=0 && rect.height()!=0)
			{
				AreaType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
				changeVisibleArea_(area, true, true);
			}
		}
		else if (action_mode_==AM_MEASURE)
		{
			measurement_start_.clear();
			update_(__PRETTY_FUNCTION__);
		}

		//do the normal stuff
		SpectrumCanvas::keyReleaseEvent(e);
	}

	void Spectrum2DCanvas::mouseDoubleClickEvent(QMouseEvent* e)
	{
		if (getCurrentLayer_().type==LayerData::DT_FEATURE)
		{
			if (selected_peak_.isValid()) //edit existing feature
			{
				FeatureEditDialog dialog(this);
				dialog.setFeature(getCurrentLayer_().features[selected_peak_.peak]);
				if (dialog.exec())
				{
					getCurrentLayer_().features[selected_peak_.peak] = dialog.getFeature();
					update_buffer_ = true;
					update_(__PRETTY_FUNCTION__);

					modificationStatus_(activeLayerIndex(), true);
				}
			}
			else //create new feature
			{
				Feature tmp;
				tmp.setRT(widgetToData_(e->pos())[1]);
				tmp.setMZ(widgetToData_(e->pos())[0]);
				FeatureEditDialog dialog(this);
				dialog.setFeature(tmp);
				if (dialog.exec())
				{
					tmp = dialog.getFeature();
					getCurrentLayer_().features.push_back(tmp);
					update_buffer_ = true;
					update_(__PRETTY_FUNCTION__);

					modificationStatus_(activeLayerIndex(), true);
				}
			}
		}
	}

	void Spectrum2DCanvas::mergeIntoLayer(Size i, FeatureMapType& map)
	{
		OPENMS_PRECONDITION(i < layers_.size(), "Spectrum2DCanvas::mergeIntoLayer(i, map) index overflow");
		OPENMS_PRECONDITION(layers_[i].type==LayerData::DT_FEATURE, "Spectrum2DCanvas::mergeIntoLayer(i, map) non-feature layer selected");
		//reserve enough space
		layers_[i].features.reserve(layers_[i].features.size()+map.size());
		//add features
		for (Size j=0; j<map.size(); ++j)
		{
			layers_[i].features.push_back(map[j]);
		}
		//update the layer and overall ranges (if necessary)
		RangeManager<2>::PositionType min_pos_old = layers_[i].features.getMin();
		RangeManager<2>::PositionType max_pos_old = layers_[i].features.getMax();
		DoubleReal min_int_old = layers_[i].features.getMinInt();
		DoubleReal max_int_old = layers_[i].features.getMaxInt();
		layers_[i].features.updateRanges();
		if(min_pos_old>layers_[i].features.getMin() || max_pos_old<layers_[i].features.getMax())
		{
			recalculateRanges_(0,1,2);
			resetZoom(true);
		}
		if(min_int_old>layers_[i].features.getMinInt() || max_int_old<layers_[i].features.getMaxInt())
		{
			intensityModeChange_();
		}
	}

	void Spectrum2DCanvas::mergeIntoLayer(Size i, ConsensusMapType& map)
	{
		OPENMS_PRECONDITION(i < layers_.size(), "Spectrum2DCanvas::mergeIntoLayer(i, map) index overflow");
		OPENMS_PRECONDITION(layers_[i].type==LayerData::DT_CONSENSUS, "Spectrum2DCanvas::mergeIntoLayer(i, map) non-consensus-feature layer selected");
		//reserve enough space
		layers_[i].consensus.reserve(layers_[i].features.size()+map.size());
		//add features
		for (Size j=0; j<map.size(); ++j)
		{
			layers_[i].consensus.push_back(map[j]);
		}
		//update the layer and overall ranges (if necessary)
		RangeManager<2>::PositionType min_pos_old = layers_[i].consensus.getMin();
		RangeManager<2>::PositionType max_pos_old = layers_[i].consensus.getMax();
		DoubleReal min_int_old = layers_[i].consensus.getMinInt();
		DoubleReal max_int_old = layers_[i].consensus.getMaxInt();
		layers_[i].consensus.updateRanges();
		if(min_pos_old>layers_[i].consensus.getMin() || max_pos_old<layers_[i].consensus.getMax())
		{
			recalculateRanges_(0,1,2);
			resetZoom(true);
		}
		if(min_int_old>layers_[i].consensus.getMinInt() || max_int_old<layers_[i].consensus.getMaxInt())
		{
			intensityModeChange_();
		}
	}

} //namespace OpenMS

