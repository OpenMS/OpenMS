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


// Qt
#include <QtGui/QResizeEvent>
#include <QtGui/QMouseEvent>
#include <QtGui/QPaintEvent>
#include <QtGui/QPainter>

// STL
#include <iostream>

// OpenMS
#include <OpenMS/VISUAL/HistogramWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	
	HistogramWidget::HistogramWidget(const Histogram<UnsignedInt,float>& distribution, QWidget* parent)
	  : QWidget(parent),
		dist_(distribution),
		show_splitters_(false),
		moving_splitter_(0),
		margin_(30),
		buffer_(),
		scaling_factor_(100)
	{
     //use log scale for int 	 
     dist_.applyLogTransformation(100);
     
		left_splitter_ =  dist_.min();
		right_splitter_ = dist_.max();
		setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
		setMinimumSize(600,450);
		bottom_axis_ = new AxisWidget(AxisWidget::BOTTOM,"x",this);
		//bottom_axis_->setPaletteBackgroundColor(Qt::yellow); //for debugging:
		bottom_axis_->setMargin(margin_);
		bottom_axis_->setTickLevel(2);
		bottom_axis_->setAxisBounds(dist_.min(),dist_.max());
	}
	
	HistogramWidget::~HistogramWidget()
	{
		delete(bottom_axis_);
	}
	
	float HistogramWidget::getLeftSplitter()
	{
		return left_splitter_;	
	}
	
	float HistogramWidget::getRightSplitter()
	{
		return right_splitter_;	
	}
	
	void HistogramWidget::showSplitters(bool on)
	{
		show_splitters_=on;	
	}
	
	void HistogramWidget::setRightSplitter(float pos)
	{
		right_splitter_=min(dist_.max(),pos);
	}
	
	void HistogramWidget::setLeftSplitter(float pos)
	{
		left_splitter_=max(dist_.min(),pos);
	}
	
	void HistogramWidget::setLegend(const string& legend)
	{
		bottom_axis_->setLegend(legend);
	}
	
	
	void HistogramWidget::mousePressEvent( QMouseEvent *e)
	{
		if (show_splitters_ && e->button()==Qt::LeftButton)
		{
			//left
			SignedInt p = margin_ + UnsignedInt(((left_splitter_-dist_.min())/(dist_.max()-dist_.min()))*(width()-2*margin_));
			//cout << "Mouse: " << e->x() << " p: " << p << " splitter: " << left_splitter_ << endl;
			if (e->x()>=p && e->x()<=p+5)
			{
				moving_splitter_=1;
			}
			
			//right
			p = margin_ + UnsignedInt(((right_splitter_-dist_.min())/(dist_.max()-dist_.min()))*(width()-2*margin_));
			if (e->x()<=p && e->x()>=p-5)
			{
				moving_splitter_=2;
			}
		}
		else
		{
			e->ignore();
		}
	}
	
	void HistogramWidget::mouseMoveEvent( QMouseEvent *e)
	{
		if (show_splitters_ && (e->buttons() & Qt::LeftButton))
		{
			//left
			if (moving_splitter_==1)
			{
				left_splitter_ = float(SignedInt(e->x())-SignedInt(margin_))/(width()-2*margin_)*(dist_.max()-dist_.min())+dist_.min();
				//upper bound
				if (left_splitter_>right_splitter_-(dist_.max()-dist_.min())/50.0)
				{
					left_splitter_ = right_splitter_-(dist_.max()-dist_.min())/50.0;
				}
				//lower bound
				if (left_splitter_<dist_.min()) 
				{
					left_splitter_=dist_.min();
				}
				update();
			}
			
			//right
			if (moving_splitter_==2)
			{
				
				right_splitter_ = float(SignedInt(e->x())-SignedInt(margin_))/(width()-2*margin_+2)*(dist_.max()-dist_.min())+dist_.min();
				//upper bound
				if (right_splitter_<left_splitter_+(dist_.max()-dist_.min())/50.0)
				{
					right_splitter_ = left_splitter_+(dist_.max()-dist_.min())/50.0;
				}
				//lower bound
				if (right_splitter_>dist_.max()) 
				{
					right_splitter_=dist_.max();
				}
				update();
			}
		}
		else
		{
			e->ignore();
		}
	}
	
	void HistogramWidget::mouseReleaseEvent( QMouseEvent *e)
	{
		if (show_splitters_)
		{
			moving_splitter_ = 0;
		}
		else
		{
			e->ignore();
		}
	}
	
	void HistogramWidget::paintEvent( QPaintEvent * /*e*/)
	{
		QPainter painter2(this);
		painter2.drawPixmap(margin_, 0, buffer_);
		
		//draw splitters
		if (show_splitters_)
		{
			QPainter painter(this);
			painter.setPen(Qt::green);
			
			//cout << "Left splitter: " << left_splitter_<< " dist: " << dist_.min() << endl;
			//cout << "Right splitter: " << right_splitter_<< " dist: " << dist_.max() << endl;
	
			//left
			UnsignedInt p =  UnsignedInt(((left_splitter_-dist_.min())/(dist_.max()-dist_.min()))*(width()-2*margin_))+margin_;
			//cout << "Left splitter position: " << p << endl;
			painter.drawLine(p,margin_-8,p,height()-bottom_axis_->height());
			painter.drawLine(p,margin_-8,p+5,margin_-8);
			painter.drawLine(p+5,margin_-8,p,margin_-3);
	
			//right
			p = UnsignedInt(((right_splitter_-dist_.min())/(dist_.max()-dist_.min()))*(width()-2*margin_))+margin_;
			painter.drawLine(p,margin_-8,p,height()-bottom_axis_->height());
			painter.drawLine(p,margin_-8,p-5,margin_-8);
			painter.drawLine(p-5,margin_-8,p,margin_-3);
		}
	}
	
	void HistogramWidget::resizeEvent( QResizeEvent * /*e*/)
	{
		buffer_ = QPixmap(width()-margin_,height()-bottom_axis_->height());
		bottom_axis_->setGeometry(margin_,height()-bottom_axis_->height(),width()-margin_,bottom_axis_->height());
		invalidate_();
	}
	
	void HistogramWidget::invalidate_()
	{
		QPainter painter(&buffer_);
		buffer_.fill(palette().window().color());
		UnsignedInt w = buffer_.width();
		UnsignedInt h = buffer_.height();
	
		//draw distribution	
		QPen pen;
		pen.setWidth(UnsignedInt(rint(float(w)/(dist_.size()*2)))); //reconfigure pen width
		pen.setColor(QColor(100,125,175));
		painter.setPen(pen);
		
		for (UnsignedInt i=0; i<dist_.size();++i)
		{
			UnsignedInt p = UnsignedInt((float(i)/(dist_.size()-1))*(w-margin_));
			UnsignedInt top = UnsignedInt(((dist_[i] / scaling_factor_)/(dist_.maxValue() / scaling_factor_))*(h-margin_));
			painter.drawLine(p+1,h,p+1,h-top);
		}
	
		//calculate total intensity
		float total_sum=0;
		for (UnsignedInt i=0; i<dist_.size();++i)
		{
			total_sum += (dist_[i] / scaling_factor_);
		}	
	
		// draw part of total intensity
		painter.setPen(Qt::red);
		QPoint last_point(1,h);
		QPoint point;
		float int_sum=0;
		for (UnsignedInt i=0; i<dist_.size();++i)
		{
			int_sum += (dist_[i] / scaling_factor_);
			point.setX(UnsignedInt((float(i)/(dist_.size()-1))*(w-margin_)));
			point.setY(UnsignedInt((1-(int_sum / total_sum))*(h-margin_)+margin_));
			painter.drawLine(last_point,point);
			last_point=point;
		}	
		//draw coord system	(on top distribution because of pen width)
		painter.setPen(Qt::black);
	
		painter.drawLine(0,margin_,0,h-1);
		painter.drawLine(0,h-1,w-margin_+1,h-1);
		
		update();
	}

} //namespace OpenMS
