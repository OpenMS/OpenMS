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
		action_mode_ = AM_TRANSLATE;
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
	
	SignedInt Spectrum3DCanvas::finishAdding(float low_intensity_cutoff)
	{
		if (layers_.back().type!=LayerData::DT_PEAK)
		{
			return -1;
		}
		
		current_layer_ = getLayerCount()-1;
		currentPeakData_().sortSpectra(true);
		currentPeakData_().updateRanges(1);	
		recalculateRanges_(1,0,2);
		//values for datareduction
		sum_of_peaks_ = getCurrentPeakData().getSize();
		area_ = (getCurrentPeakData().getMaxRT()-getCurrentPeakData().getMinRT())*(getCurrentPeakData().getMaxMZ()-getCurrentPeakData().getMinMZ());
		peaks_per_rt_ = (int)floor(sum_of_peaks_/getCurrentPeakData().size());
	
	 	if(getPrefAsString("Preferences:3D:Reduction:Mode")!="Reduction OFF")
	 	{
			makeReducedDataSet();
	 	}
	
		visible_area_.assign(overall_data_range_);
		
		layers_.back().min_int = low_intensity_cutoff;
		layers_.back().max_int = overall_data_range_.max_[2];
		emit layerActivated(this);
		openglwidget()->recalculateDotGradient_();
		update_(__PRETTY_FUNCTION__);
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
		if(sum_of_peaks_ < getPrefAsInt("Preferences:3D:DisplayedPeaks")||getPrefAsString("Preferences:3D:Reduction:Mode")=="Reduction OFF")
		{	
	 		show_reduced_ = false;
	 		recalculateRanges_(1,0,2);
		}
		else 
		{
			Param reduction_param;
			show_reduced_ = true;
			if(getPrefAsString("Preferences:3D:Reduction:Mode")=="MaxReduction")
			{	
				int reduction;
				if(zoom_stack_.empty())
				{
					reduction = sum_of_peaks_/getPrefAsInt("Preferences:3D:DisplayedPeaks");
				}
				else
				{
					double new_area = (visible_area_.max_[0]-visible_area_.min_[0])*(visible_area_.max_[1]-visible_area_.min_[1]);
					int needed_peak_number = (int)(getPrefAsInt("Preferences:3D:DisplayedPeaks")* area_ / new_area);
					if(needed_peak_number<sum_of_peaks_)
					{
						reduction = sum_of_peaks_/needed_peak_number;
						if (datareducer_==0 || datareducer_->getName()!="MaxReducer")
						{
							datareducer_ = Factory<DataReducer>::create("MaxReducer");
						}
					}
					else
					{
						show_reduced_= false;
					}
				}
				reduction_param.setValue("Peaksperstep", reduction);
			}
			else if(getPrefAsString("Preferences:3D:Reduction:Mode")=="SumReduction")
			{	
				double reduction;
				if(zoom_stack_.empty())
				{
					reduction = (double)sum_of_peaks_/(	(double)peaks_per_rt_*(double)getPrefAsInt("Preferences:3D:DisplayedPeaks"));
				}
				else
				{
					double new_area = (visible_area_.max_[0]-visible_area_.min_[0])*(visible_area_.max_[1]-visible_area_.min_[1]);
					int needed_peak_number = (int)(getPrefAsInt("Preferences:3D:DisplayedPeaks")* area_ / new_area);
	 				if(needed_peak_number<sum_of_peaks_)
					{
						reduction = (double)sum_of_peaks_/(	(double)peaks_per_rt_* needed_peak_number);
						if (datareducer_==0 || datareducer_->getName()!="SumReducer")
						{
							datareducer_ = Factory<DataReducer>::create("SumReducer");
						}
					}
					else
					{
						show_reduced_= false;
					}
				}
				reduction_param.setValue("Rangeperstep",reduction);
			}
			if(show_reduced_)
			{
				for(UnsignedInt i = 0; i<layers_.size();i++)
				{
					/// @todo added fix for segfault for release (Cornelia)
					if (datareducer_ != 0)
					{
						datareducer_->applyReduction(getLayer(i).peaks,getLayer_(i).reduced);
					}
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
	
	void Spectrum3DCanvas::actionModeChange_()
	{
		switch(action_mode_)
		{
			case AM_TRANSLATE:
				openglwidget()->setAngels(220,220,0);
				openglwidget()->setZoomFactor(1.5,false);
				update_buffer_ = true;
				update_(__PRETTY_FUNCTION__);
				break;
			case AM_ZOOM:
				openglwidget()->setAngels(1440,0,0);
				openglwidget()->resetTranslation();
				openglwidget()->setZoomFactor(1.25,false);
				update_buffer_ = true;
				update_(__PRETTY_FUNCTION__);
				break;
			default:
				throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				break;
			}
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
		  
		SignedInt Spectrum3DCanvas::getDataMode()
		{
			if(prefs_.getValue("Preferences:3D:Data:Mode").isEmpty())
			{
				return 0;
			}
			return SignedInt(prefs_.getValue("Preferences:3D:Data:Mode"));
		}
		
		void Spectrum3DCanvas::setDataMode()
		{
			makeReducedDataSet();
			if(zoom_stack_.empty())
			{
				resetZoom();
			}
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
		
		SignedInt Spectrum3DCanvas::getDotMode()
		{
			if (prefs_.getValue("Preferences:3D:Dot:Mode").isEmpty())
			{
				return 0;
			}
			return SignedInt(prefs_.getValue("Preferences:3D:Dot:Mode"));
		}
		
		void Spectrum3DCanvas::setDotGradient(const std::string& gradient)
		{
			openglcanvas_->setDotGradient(gradient);
		}
		
		SignedInt Spectrum3DCanvas::getShadeMode()
		{
			if(prefs_.getValue("Preferences:3D:Shade:Mode").isEmpty())
			{
				return 0;
			}
			return SignedInt(prefs_.getValue("Preferences:3D:Shade:Mode"));
		}
		
		UnsignedInt Spectrum3DCanvas::getDotInterpolationSteps()
		{
			if(prefs_.getValue("Preferences:3D:InterpolationSteps").isEmpty())
			{
				return 0;
			}
			return UnsignedInt(prefs_.getValue("Preferences:3D:InterpolationSteps"));
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
				openglwidget()->glInit(); 
			}
			openglwidget()->resizeGL(width(),height());
			openglwidget()->glDraw(); 
		}

}//namspace

