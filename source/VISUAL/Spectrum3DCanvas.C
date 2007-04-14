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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum3DCanvasPDP.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FILTERING/DATAREDUCTION/DataReducer.h>

#include <QtGui/QResizeEvent>

using namespace std;

namespace OpenMS
{
	using namespace Internal;
	
	Spectrum3DCanvas::Spectrum3DCanvas(QWidget* parent)
		: SpectrumCanvas(parent)
	{  
		setFocusPolicy(Qt::TabFocus);
		openglcanvas_= new Spectrum3DOpenGLCanvas(this, *this);
		action_mode_ = AM_ZOOM;
		legend_shown_ = true;
	}
		
	Spectrum3DCanvas::~Spectrum3DCanvas()
	{
	
	}
	
	void Spectrum3DCanvas::resizeEvent(QResizeEvent *e)
	{
		openglcanvas_ ->resize(e->size().width(),e->size().height());
	}
	
	void Spectrum3DCanvas::showLegend(bool show)
	{
		legend_shown_ = show;
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	
	Int Spectrum3DCanvas::finishAdding(float low_intensity_cutoff)
	{
		if (layers_.back().type!=LayerData::DT_PEAK)
		{
			return -1;
		}
		
		current_layer_ = getLayerCount()-1;
		currentPeakData_().sortSpectra(true);
		currentPeakData_().updateRanges(1);	
		recalculateRanges_(1,0,2);
		area_ = (getCurrentPeakData().getMaxRT()-getCurrentPeakData().getMinRT())*(getCurrentPeakData().getMaxMZ()-getCurrentPeakData().getMinMZ());
	
	 	if(getPrefAsString("Preferences:3D:ReductionMode")!="Off")
	 	{
			makeReducedDataSet();
	 	}
	
		visible_area_.assign(overall_data_range_);
		
		layers_.back().min_int = low_intensity_cutoff;
		layers_.back().max_int = overall_data_range_.max_[2];
		emit layerActivated(this);
		openglwidget()->recalculateDotGradient_();
		//	update_(__PRETTY_FUNCTION__);
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
		return current_layer_;
	}
	
	void Spectrum3DCanvas::changeVisibleArea_(const AreaType& new_area, bool add_to_stack)
	{
		if (new_area==visible_area_)
		{
			return;
		}
		//store old zoom state
		if (add_to_stack)
		{
			zoom_stack_.push(visible_area_);
		}
		visible_area_ = new_area;
		
		updateScrollbars_();
		
		emit visibleAreaChanged(new_area);
		makeReducedDataSet();
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	void Spectrum3DCanvas::makeReducedDataSet()
	{
		if(getPrefAsString("Preferences:3D:ReductionMode")=="Off")
		{	
	 		show_reduced_ = false;
	 		recalculateRanges_(1,0,2);
			return;
		}
		else if(getCurrentLayer().peaks.getSize() < (UInt)(getPrefAsInt("Preferences:3D:DisplayedPeaks")))
		{
			show_reduced_ = false;
			recalculateRanges_(1,0,2);
			return;
		}
		else
		{
			Param reduction_param;
			show_reduced_ = true;
			if(getPrefAsString("Preferences:3D:ReductionMode")=="Max reduction")
			{	
				int reduction;
				if(zoom_stack_.empty())
				{
					reduction = getCurrentLayer().peaks.getSize()/(UInt)(getPrefAsInt("Preferences:3D:DisplayedPeaks"));
				}
				else
				{
					double new_area = (visible_area_.max_[0]-visible_area_.min_[0])*(visible_area_.max_[1]-visible_area_.min_[1]);
					UInt needed_peak_number = (UInt)(getPrefAsInt("Preferences:3D:DisplayedPeaks")* area_ / new_area);
					if(needed_peak_number<getCurrentLayer().peaks.getSize())
					{
						reduction = getCurrentLayer().peaks.getSize()/needed_peak_number;
					}
					else
					{
						show_reduced_= false;
						return;
					}
				}
				if (datareducer_!= 0 && datareducer_->getName()!="MaxReducer")
				{
					delete datareducer_;
					datareducer_ = NULL;
				}
				if(datareducer_==0)
				{
					datareducer_ = Factory<DataReducer>::create("MaxReducer");
					
				}
				reduction_param.setValue("Peaksperstep", reduction);
				}
			else if(getPrefAsString("Preferences:3D:ReductionMode")=="Sum reduction")
			{	
				int peaks_per_rt = (int)floor(getCurrentLayer().peaks.getSize()/getCurrentLayer().peaks.size());
				double reduction;
				if(zoom_stack_.empty())
				{
					reduction = (double)getCurrentLayer().peaks.getSize()/(	(double)peaks_per_rt*(double)getPrefAsInt("Preferences:3D:DisplayedPeaks"));
				}
				else
				{
					double new_area = (visible_area_.max_[0]-visible_area_.min_[0])*(visible_area_.max_[1]-visible_area_.min_[1]);
					UInt needed_peak_number = (UInt)(getPrefAsInt("Preferences:3D:DisplayedPeaks")* area_ / new_area);
	 				if(needed_peak_number<getCurrentLayer().peaks.getSize())
					{
						reduction = (double)getCurrentLayer().peaks.getSize()/(	(double)peaks_per_rt* needed_peak_number);
					}
						else
					{
						show_reduced_= false;
						return;
					}
				}
				if (datareducer_!=0 && datareducer_->getName()!="SumReducer")
				{
					delete datareducer_;
					datareducer_ = NULL;
				}
				if (datareducer_==0)
				{
					datareducer_ = Factory<DataReducer>::create("SumReducer");
					}
				reduction_param.setValue("Rangeperstep",reduction);
			}
			
			if(show_reduced_)
			{
				for(UInt i = 0; i<layers_.size();i++)
				{
						datareducer_->applyReduction(getLayer(i).peaks,getLayer_(i).reduced);
				}
				recalculateRanges_(1,0,2);
			}
		}
	}
	
	void Spectrum3DCanvas::setMainPreferences(const Param& prefs)
	{
		SpectrumCanvas::setMainPreferences(prefs);
		openglwidget()->gradient_.fromString(getPrefAsString("Preferences:3D:Dot:Gradient"));
	}
	
	PreferencesDialogPage * Spectrum3DCanvas::createPreferences(QWidget* parent)
	{
		return new Spectrum3DCanvasPDP(this,parent);
	}

	void Spectrum3DCanvas::activateLayer(int layer_index)
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()) || layer_index==int(current_layer_))
		{
			return ;
		}
		current_layer_ = layer_index;
		emit layerActivated(this);
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum3DCanvas::repaintAll()
	{
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum3DCanvas::intensityModeChange_()
	{
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum3DCanvas::removeLayer(int layer_index)
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()))
		{
			return;
		}
		layers_.erase(layers_.begin()+layer_index);
		
		//update current layer
		if (current_layer_!=0 && current_layer_ >= getLayerCount())
		{
		current_layer_ = getLayerCount()-1;
		}
		
		recalculateRanges_(1,0,2);
		visible_area_.assign(overall_data_range_);
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	
	Spectrum3DOpenGLCanvas* Spectrum3DCanvas::openglwidget()
	{
		return static_cast<Spectrum3DOpenGLCanvas*>(openglcanvas_);
	}
	
	
	////preferences////////////////////
	
	void Spectrum3DCanvas::setDotGradient(const std::string& gradient)
	{
		openglcanvas_->setDotGradient(gradient);
	}
	
	Int Spectrum3DCanvas::getShadeMode()
	{
		if(prefs_.getValue("Preferences:3D:Dot:ShadeMode").isEmpty())
		{
			return 0;
		}
		return Int(prefs_.getValue("Preferences:3D:Dot:ShadeMode"));
	}

	void Spectrum3DCanvas::update_(const char*
#ifdef DEBUG_UPDATE_
			caller_name)
	{
		cout << "Spectrum3DCanvas::update_ from '" << caller_name << "'" << endl;
#else
		)
	{
#endif
		if(update_buffer_)
		{
			update_buffer_ = false;
			if(intensity_mode_ == SpectrumCanvas::IM_SNAP)
			{
				openglwidget()->updateIntensityScale();
			}
			openglwidget()->recalculateDotGradient_();
			openglwidget()->initializeGL(); 
		}
		openglwidget()->resizeGL(width(),height());	
		openglwidget()->glDraw(); 
		
	}

}//namspace

