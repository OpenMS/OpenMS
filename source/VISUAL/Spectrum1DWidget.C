// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: Spectrum1DWidget.C,v 1.36 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// Qt
#include <qaction.h>

// OpenMS
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DWidgetPDP.h>


using namespace std;

namespace OpenMS
{
	using namespace Internal;
	
	Spectrum1DWidget::Spectrum1DWidget(QWidget* parent, const char* name, WFlags f)
		: SpectrumWidget(parent, name, f),
		label_mode_(Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT),
		log_switched_(false)
	{
		// TODO: TICWidget overwrites 1DCanvas (canvas_)
		//set the label mode for the axes  - side effect
		setCanvas(new Spectrum1DCanvas(this, "Spectrum1DCanvas"));
		connect(canvas(), SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)),
		        this, SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)));
		connect(canvas(), SIGNAL(sendCursorStatus(double,double,double)),
		        this, SIGNAL(sendCursorStatus(double,double,double)));
		setLabelMode_(label_mode_);
		x_axis_->setLegend("m/z");
		y_axis_->setLegend("Intensity");
		addClient(canvas(),"Canvas",true);
	
		setMouseTracking(true);
		connect(x_axis_, SIGNAL(sendAxisMouseMovement()),this, SLOT(clearHighlighting()));
		connect(y_axis_, SIGNAL(sendAxisMouseMovement()),this, SLOT(clearHighlighting()));
	}
	
	Spectrum1DCanvas* Spectrum1DWidget::canvas() const
	{
		return static_cast<Spectrum1DCanvas*>(canvas_);
	}
	
	void Spectrum1DWidget::setMainPreferences(const Param& prefs)
	{
		SpectrumWidget::setMainPreferences(prefs);
		
		x_axis_->showLegend(getPrefAsInt("Preferences:1D:X:Legend")==1);
		y_axis_->showLegend(getPrefAsInt("Preferences:1D:Y:Legend")==1);
		
		if (getPrefAsInt("Preferences:1D:Intensity:Logarithmic")==1) setIntensityModificationLog();
		if (getPrefAsInt("Preferences:1D:Intensity:Relative")==0) intensityAxisAbsolute();
	}
	
	
	void Spectrum1DWidget::mouseMoveEvent( QMouseEvent* /*e*/)
	{
		clearHighlighting();
	}
	
	void Spectrum1DWidget::clearHighlighting()
	{
		canvas()->clearHighlighting();
	}
	
	void Spectrum1DWidget::setVisibleArea(double x1, double x2)
	{
		canvas()->setVisibleArea(x1, x2);
	}
	
	void Spectrum1DWidget::switchAxis(bool swapped_axes) // TODO: fix name and typo. should be setSwappedAxes()
	{
		switch (label_mode_)
		{
			case Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT:
				setLabelMode_(Spectrum1DCanvas::LM_XPERCENT_YABSOLUTE);
			break;
			case Spectrum1DCanvas::LM_XPERCENT_YABSOLUTE:
				setLabelMode_(Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT);
			break;
			default:  // in the other 2 cases nothing needs to be done
			break;
		}
	
		SpectrumWidget::switchAxis(swapped_axes);
	}
	
	void Spectrum1DWidget::recalculateAxes()
	{
		
		const SpectrumCanvas::AreaType& visible_area = canvas()->visible_area_;
		const SpectrumCanvas::AreaType& data_area = canvas()->getDataArea_();
		
		//cout << "Spectrum1DWidget::recalculateAxes() IN(data): x: " << visible_area.minX() << " " <<visible_area.maxX() << "   y: " <<visible_area.minY() << " " <<visible_area.maxY() << endl;
		//cout << "Spectrum1DWidget::recalculateAxes() IN(visible): x: " << data_area.minX() << " " <<data_area.maxX() << "   y: " <<data_area.minY() << " " <<data_area.maxY() << endl;
	
		const MappingInfo& mapping_info = *canvas()->getMappingInfo();
	
		//calculate margins around data_area
		double mx = 0.002*(data_area.maxX() - data_area.minX());
		double my = 0.002*(data_area.maxY() - data_area.minY());
	
		//fixed bounds for relative display
		double rel_lo = 100.0/data_area.maxY();
		if (canvas()->getIntensityModification()==SpectrumCanvas::IM_LOG)
			rel_lo = 1.0-(0.002*100.0);
		double rel_hi = 100+(0.002*100.0);
	
		// recalculate gridlines
		double lx,hx,ly,hy;
	
		switch (label_mode_)
		{
			case Spectrum1DCanvas::LM_XABSOLUTE_YABSOLUTE:
			if (mapping_info.isMzToXAxis())
			{
				lx = visible_area.minX();
				hx = visible_area.maxX();
				ly = visible_area.minY();
				hy = visible_area.maxY();
			}
			else
			{
				lx = visible_area.minY();
				hx = visible_area.maxY();
				ly = visible_area.minX();
				hy = visible_area.maxX();
			}
			break;
	
			case Spectrum1DCanvas::LM_XPERCENT_YABSOLUTE:
			if (mapping_info.isMzToXAxis())
			{
				lx = intervalTransformation(visible_area.minX(), data_area.minX()-mx, data_area.maxX()+mx, rel_lo, rel_hi);
				hx = intervalTransformation(visible_area.maxX(), data_area.minX()-mx, data_area.maxX()+mx, rel_lo, rel_hi);
				ly = visible_area.minY();
				hy = visible_area.maxY();
			} else
			{
				lx = intervalTransformation(visible_area.minY(), data_area.minX()-mx, data_area.maxX()+mx, rel_lo, rel_hi);
				hx = intervalTransformation(visible_area.maxY(), data_area.minX()-mx, data_area.maxX()+mx, rel_lo, rel_hi);
				ly = visible_area.minX();
				hy = visible_area.maxX();
			}
			break;
	
			case Spectrum1DCanvas::LM_XPERCENT_YPERCENT:
			if (mapping_info.isMzToXAxis())
			{
				lx = intervalTransformation(visible_area.minX(), data_area.minX()-mx, data_area.maxX()+mx, rel_lo, rel_hi);
				hx = intervalTransformation(visible_area.maxX(), data_area.minX()-mx, data_area.maxX()+mx, rel_lo, rel_hi);
				ly = intervalTransformation(visible_area.minY(), data_area.minY()-my, data_area.maxY()+my, rel_lo, rel_hi);
				hy = intervalTransformation(visible_area.maxY(), data_area.minY()-my, data_area.maxY()+my, rel_lo, rel_hi);
			} else
			{
				lx = intervalTransformation(visible_area.minY(), data_area.minY()-my, data_area.maxY()+my, rel_lo, rel_hi);
				hx = intervalTransformation(visible_area.maxY(), data_area.minY()-my, data_area.maxY()+my, rel_lo, rel_hi);
				ly = intervalTransformation(visible_area.minX(), data_area.minX()-mx, data_area.maxX()+mx, rel_lo, rel_hi);
				hy = intervalTransformation(visible_area.maxX(), data_area.minX()-mx, data_area.maxX()+mx, rel_lo, rel_hi);
			}
			break;
	
			default:	// LM_XABSOLUTE_YPERCENT
			if (mapping_info.isMzToXAxis())
			{
				lx = visible_area.minX();
				hx = visible_area.maxX();
				ly = intervalTransformation(visible_area.minY(), data_area.minY()-my, data_area.maxY()+my, rel_lo, rel_hi);
				hy = intervalTransformation(visible_area.maxY(), data_area.minY()-my, data_area.maxY()+my, rel_lo, rel_hi);
			} else
			{
				lx = visible_area.minY();
				hx = visible_area.maxY();
				ly = intervalTransformation(visible_area.minX(), data_area.minX()-mx, data_area.maxX()+mx, rel_lo, rel_hi);
				hy = intervalTransformation(visible_area.maxX(), data_area.minX()-mx, data_area.maxX()+mx, rel_lo, rel_hi);
			}
			break;
		}
		//cout << "Spectrum1DWidget::recalculateAxes() OUT: x: " << lx << " " <<hx << "   y: " <<ly << " " <<hy << endl;
		
		// Undo logarithmic scale, setAxis deals with plain min and max
		bool is_intensity_axis_percent = !canvas()->isAbsoluteIntensity();
		if ((canvas()->getIntensityModification() == SpectrumCanvas::IM_NONE && log_switched_)
				|| (canvas()->getIntensityModification() == SpectrumCanvas::IM_LOG && !log_switched_))
		{
			if (mapping_info.isMzToXAxis()) // y = intensity
			{
				if (canvas()->getSnapToMax()) // Adjust axis values in snap-to-max-intensity-mode
					y_axis_->setAxisBounds(log2linear(ly/canvas()->getSnapFactor(), is_intensity_axis_percent, old_max_intensity_),
																	 log2linear(hy/canvas()->getSnapFactor(), is_intensity_axis_percent, old_max_intensity_));
				else
					y_axis_->setAxisBounds(log2linear(ly, is_intensity_axis_percent, old_max_intensity_),
																	 log2linear(hy, is_intensity_axis_percent, old_max_intensity_));
			}
			else // x = intensity
			{
				if (canvas()->getSnapToMax()) // Adjust axis values in snap-to-max-intensity-mode
					x_axis_->setAxisBounds(log2linear(lx/canvas()->getSnapFactor(), is_intensity_axis_percent, old_max_intensity_),
																	 log2linear(hx/canvas()->getSnapFactor(), is_intensity_axis_percent, old_max_intensity_));
				else
					x_axis_->setAxisBounds(log2linear(lx, is_intensity_axis_percent, old_max_intensity_),
																	 log2linear(hx, is_intensity_axis_percent, old_max_intensity_));
			}
		}
		else  // no logarithmic scale
		{
			if (mapping_info.isMzToXAxis())  // y = intensity
			{
				x_axis_->setAxisBounds(lx, hx);
				if (canvas()->getSnapToMax()) // Adjust axis values in snap-to-max-intensity-mode
					y_axis_->setAxisBounds(ly/canvas()->getSnapFactor(), hy/canvas()->getSnapFactor());
				else
					y_axis_->setAxisBounds(ly, hy);
			}else  // x = intensity
			{
				y_axis_->setAxisBounds(ly, hy);
				if (canvas()->getSnapToMax()) // Adjust axis values in snap-to-max-intensity-mode
					x_axis_->setAxisBounds(lx/canvas()->getSnapFactor(), hx/canvas()->getSnapFactor());
				else
					x_axis_->setAxisBounds(lx, hx);
			}
		}
		log_switched_ = false;
	}
	
	void Spectrum1DWidget::setZoomFactor(double d)
	{
		canvas()->setZoomFactor(d);
	}
	
	QImage Spectrum1DWidget::getImage(UnsignedInt width, UnsignedInt height, UnsignedInt flags)
	{
		this->setUpdatesEnabled(false);
	
	 	QColor c_a = x_axis_->paletteBackgroundColor();
	 	QColor c_o = y_axis_->paletteBackgroundColor();
	 	QColor c_t = paletteBackgroundColor();
	
		y_axis_->setPaletteBackgroundColor(white);
		x_axis_->setPaletteBackgroundColor(white);
		setBackgroundColor(white);
	
		int h = this->height();
		int w = this->width();
	
		int pen=0;
		if (flags & THICKLINES)
			pen = 2+width/1024;
		else
			pen = width/1024;
	
		canvas()->setPenWidth(pen);
		y_axis_->setPenWidth(pen);
		x_axis_->setPenWidth(pen);
	
		resize(width,height);
		hide();
	
		QPixmap p = QPixmap::grabWidget(this);
		QImage image = p.convertToImage();
	
		resize(w,h);
	
	 	y_axis_->setPaletteBackgroundColor(c_o);
	 	x_axis_->setPaletteBackgroundColor(c_a);
		setBackgroundColor(c_t);
	
	 
		if (flags & GREYSCALE)
		{
			for (unsigned int y = 0; y!= height;y++)
			{
				unsigned int *p = (unsigned int*)image.scanLine(y);
				for (unsigned int x = 0;x != width;x++)
				{
					// TODO: add brightness/contrast control
					unsigned int cgrey = static_cast<int>((qRed(*p)+qGreen(*p)+qBlue(*p))/3);
					*p++ = qRgb(cgrey,cgrey,cgrey);
				}
			}
		}
	
		canvas()->setPenWidth(0);
		y_axis_->setPenWidth(0);
		x_axis_->setPenWidth(0);
	
		this->setUpdatesEnabled(true);
		y_axis_->repaint();
	 	x_axis_->repaint();
	
		show();
		return image;
	}
	
	void Spectrum1DWidget::setDrawMode(QAction* a)
	{
		canvas()->setDrawMode(a);
	}
	
	void Spectrum1DWidget::drawModePeaks()
	{
		canvas()->drawModePeaks();
	}
	
	void Spectrum1DWidget::drawModeLines()
	{
		canvas()->drawModeLines();
	}
	
	void Spectrum1DWidget::intensityModificationChange_()
	{
		log_switched_ = true;
		//cout << "Spectrum1DWidget::intensityModificationChange_()"<<endl;
		setLabelMode_(label_mode_);  // update label mode of axes (intensity_modification_ is used there)
		//cout << "Spectrum1DWidget::intensityModificationChange_() NACH LABEL MODE"<<endl;
	 	canvas()->intensityModificationChange_();
		if (getSnapToMax())	setSnapToMax(true);
		recalculateAxes();
	}
	
	Histogram<UnsignedInt,float> Spectrum1DWidget::createIntensityDistribution_()
	{
		Histogram<UnsignedInt,float> tmp(canvas()->getCurrentMinIntensity(),canvas()->getCurrentMaxIntensity(),(canvas()->getCurrentMaxIntensity() - canvas()->getCurrentMinIntensity())/500);
	
		for (Spectrum1DCanvas::ExperimentType::SpectrumType::ConstIterator it = canvas()->currentDataSet()[0].begin(); it != canvas()->currentDataSet()[0].end(); ++it)
		{
			tmp.inc(it->getIntensity());
		}
		return tmp;
	}
	
	void Spectrum1DWidget::legendModificationChange_()
	{
		y_axis_->showLegend(show_legend_);
		x_axis_->showLegend(show_legend_);
		update();
	}
	
	void Spectrum1DWidget::intensityAxisAbsolute()
	{
		bool is_log = isLogIntensity();
		if (is_log) setIntensityModificationNone();
		setIntensityAxisAbsolute_();
	 	canvas()->intensityAxisAbsolute();
		if (is_log) setIntensityModificationLog();
		if (getSnapToMax())	setSnapToMax(true);
	}
	
	void Spectrum1DWidget::intensityAxisRelative()
	{
		bool is_log = isLogIntensity();
		if (is_log) setIntensityModificationNone();
		setIntensityAxisRelative_();
		canvas()->intensityAxisRelative();
		if (is_log) setIntensityModificationLog();
		if (getSnapToMax())	setSnapToMax(true);
	}
	
	bool Spectrum1DWidget::isAbsoluteIntensity() const
	{
		return canvas()->isAbsoluteIntensity();
	}
	
	// destructor
	Spectrum1DWidget::~Spectrum1DWidget()
	{
		
	}
	
	//  label mode
	void Spectrum1DWidget::setLabelMode_(UnsignedInt label_mode)
	{
		//cout << "Spectrum1DWidget::setLabelMode_"<<endl;
		label_mode_ = label_mode;
		if (canvas()->getMappingInfo()->isMzToXAxis())
		{
			// y-axis (ordinate) is intensity axis
			x_axis_->setLogScale(false);
			y_axis_->setLogScale(isLogIntensity());
		}
		else
		{
			// x-axis is intensity axis
			y_axis_->setLogScale(false);
			x_axis_->setLogScale(isLogIntensity());
		}
		//cout << "Spectrum1DWidget::setLabelMode_ VOR recalculate"<<endl;
		recalculateAxes();
		//cout << "Spectrum1DWidget::setLabelMode_ NACH recalculate"<<endl;
	}
	
	// mapping safe wrappers start here
	void Spectrum1DWidget::setIntensityAxisAbsolute_()
	// Just a simple wrapper that simplifies using the enumeration (label_mode_) under different mappings (MappingInfo)
	{
		if (canvas()->getMappingInfo()->isMzToXAxis())
		{
			// y axis is intensity Axis
			if (label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YABSOLUTE || label_mode_ == Spectrum1DCanvas::LM_XPERCENT_YABSOLUTE) return; // intensity axis already set to absolute
			if (label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT)
			{
				setLabelMode_(Spectrum1DCanvas::LM_XABSOLUTE_YABSOLUTE);
				return;
			} else if (label_mode_ == Spectrum1DCanvas::LM_XPERCENT_YPERCENT)
			{
				setLabelMode_(Spectrum1DCanvas::LM_XPERCENT_YABSOLUTE);
				return;
			}
		} else
		{
			// x axis is intensity axis
			if (label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YABSOLUTE || label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT) return; // intensity axis already set to absolute
			if (label_mode_ == Spectrum1DCanvas::LM_XPERCENT_YABSOLUTE)
			{
				setLabelMode_(Spectrum1DCanvas::LM_XABSOLUTE_YABSOLUTE);
				return;
			} else if (label_mode_ == Spectrum1DCanvas::LM_XPERCENT_YPERCENT)
			{
				setLabelMode_(Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT);
				return;
			}
		}
	}
	
	void Spectrum1DWidget::setIntensityAxisRelative_()
	{
		if (canvas()->getMappingInfo()->isMzToXAxis())
		{
			// y axis is intensity Axis
			if (label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT || label_mode_ == Spectrum1DCanvas::LM_XPERCENT_YPERCENT) return; // intensity axis already set to relative values (%)
			if (label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YABSOLUTE)
			{
				setLabelMode_(Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT);
				return;
			} else if (label_mode_ == Spectrum1DCanvas::LM_XPERCENT_YABSOLUTE)
			{
				setLabelMode_(Spectrum1DCanvas::LM_XPERCENT_YPERCENT);
				return;
			}
		} else
		{
			// x axis is intensity axis
			if (label_mode_ == Spectrum1DCanvas::LM_XPERCENT_YABSOLUTE || label_mode_ == Spectrum1DCanvas::LM_XPERCENT_YPERCENT) return; // intensity axis already set to relative values (%)
			if (label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YABSOLUTE)
			{
				setLabelMode_(Spectrum1DCanvas::LM_XPERCENT_YABSOLUTE);
				return;
			} else if (label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT)
			{
				setLabelMode_(Spectrum1DCanvas::LM_XPERCENT_YPERCENT);
				return;
			}
		}
	}
	
	bool Spectrum1DWidget::isIntensityAxisAbsolute_() const
	{
		if (canvas()->getMappingInfo()->isMzToXAxis())
		{
			// y axis is intensity Axis
			return (label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YABSOLUTE || label_mode_ == Spectrum1DCanvas::LM_XPERCENT_YABSOLUTE);
		}
		
		// x axis is intensity axis
		return (label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YABSOLUTE || label_mode_ == Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT);
	}
	
	bool Spectrum1DWidget::getSnapToMax()
	{
		return canvas()->getSnapToMax();
	}
	
	void Spectrum1DWidget::setSnapToMax(bool b)
	{
		canvas()->setSnapToMax(b);
		recalculateAxes();
	}
	
	
	PreferencesDialogPage* Spectrum1DWidget::createPreferences(QWidget* parent)
	{
		PreferencesDialogPage* background = new Spectrum1DWidgetPDP(this, parent);
		return background;
	}
	
	void Spectrum1DWidget::minimizeToChart()
	{
		y_axis_->hide();
		x_axis_->hide();
	}

} //namespace


