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

// OpenMS includes
#include <OpenMS/VISUAL/MultiGradientSelector.h>

//qt includes
#include <QtGui/QPainter>
#include <QtGui/QColorDialog>
#include <QtGui/QPixmap>
#include <QtGui/QMouseEvent>
#include <QtGui/QKeyEvent>
#include <QtGui/QPaintEvent>

using namespace std;

namespace OpenMS
{

	MultiGradientSelector::MultiGradientSelector( QWidget * parent) 
		: QWidget(parent), 
			gradient_(), 
			margin_(5), 
			gradient_area_width_(0), 
			lever_area_height_(17), 
			selected_(-1), 
			selected_color_(Qt::white)
	{
		setMinimumSize(250,45);
		setFocusPolicy(Qt::ClickFocus);
		setToolTip( "Click the lever area to add new levers.<BR>"
								"Levers are removed with the DEL key.<BR>"
                "<NOBR>Double click a lever to change its color.</NOBR><BR>"
                "Levers can be dragged.<BR>"
                "Click the gradient to change its mode.<BR>"
		           );
	}
	
	MultiGradientSelector::~MultiGradientSelector()
	{
	
	}
	
	const MultiGradient& MultiGradientSelector::gradient() const
	{
		return gradient_;
	}
	
	MultiGradient& MultiGradientSelector::gradient()
	{
		return gradient_;
	}
	
	void MultiGradientSelector::paintEvent(QPaintEvent * /* e */ )
	{
	  static QPixmap pixmap = QPixmap(size());
		pixmap.fill(palette().window().color());	
		
		//calculate gradient area width
		if (gradient_area_width_==0)
		{
			gradient_area_width_ = width()-2*margin_-2;
		}
	
		QPainter painter(&pixmap);
		
		//gradient field outline
		painter.setPen(QColor(0,0,0));
		painter.drawRect(margin_,margin_,width()-2*margin_,height()-2*margin_-lever_area_height_);
	
		//draw gradient
		for (SignedInt i=0;i<=gradient_area_width_;++i)
		{
			painter.setPen(gradient_.interpolatedColorAt(i,0,gradient_area_width_));
			painter.drawLine(margin_+1+i,margin_+1,margin_+1+i,height()-margin_-lever_area_height_-1);		
		}
		
		//levers
		painter.setPen(QColor(0,0,0));
		for (UnsignedInt i=0;i<gradient_.size();++i)
		{
			SignedInt pos = SignedInt(float(gradient_.position(i))/100.0*gradient_area_width_+margin_+1);
			painter.drawRect(pos-4,height()-margin_-lever_area_height_+5,9,9);
			painter.drawLine(pos-4,height()-margin_-lever_area_height_+5,pos,height()-margin_-lever_area_height_);
			painter.drawLine(pos,height()-margin_-lever_area_height_,pos+4,height()-margin_-lever_area_height_+5);
			painter.fillRect(pos-3,height()-margin_-lever_area_height_+6,8,8,gradient_.color(i));
			
			//selected lever
			if (SignedInt(gradient_.position(i)) == selected_)
			{
				painter.fillRect(pos-2,height()-margin_-lever_area_height_+3,6,3,QColor(0,0,0));
				painter.fillRect(pos-1,height()-margin_-lever_area_height_+1,4,3,QColor(0,0,0));
			}
		}
		
		QPainter painter2(this);
		painter2.drawPixmap(0,0,pixmap);		
	}
	
	void MultiGradientSelector::mousePressEvent ( QMouseEvent * e )
	{
		if ( e->button() != Qt::LeftButton ) 
		{
			e->ignore();
			return;
	  } 	
		
		left_button_pressed_=true;
		
		//select lever
		for (UnsignedInt i=0;i<gradient_.size();++i)
		{
			SignedInt pos = SignedInt(float(gradient_.position(i))/100.0*gradient_area_width_+margin_+1);
			if (e->x() >= pos-3 && e->x() <= pos+4 && e->y() >= height()-margin_-lever_area_height_+8 && e->y() <= height()-margin_-lever_area_height_+15)
			{
				selected_=gradient_.position(i);
				selected_color_=gradient_.color(i);
				repaint();
				return;
			}
		}		
		
		//create new lever
		if (e->x() >= margin_ && e->x() <= width()-margin_ && e->y() >= height()-margin_-lever_area_height_ && e->y() <= height()-margin_)
		{
			SignedInt pos = SignedInt(100*(e->x()-margin_)/float(gradient_area_width_));
			gradient_.insert(pos,selected_color_);
			selected_ = pos;
			repaint();
		}
		//tmp!!
		else
		{
			if (getInterpolationMode()==MultiGradient::IM_LINEAR)
			{
				setInterpolationMode(MultiGradient::IM_STAIRS);
			}
			else
			{
				setInterpolationMode(MultiGradient::IM_LINEAR);
			}
			repaint();
		}
	}
	
		void MultiGradientSelector::mouseMoveEvent(QMouseEvent* e)
		{
			if (left_button_pressed_ && selected_!=-1)
			{
				//inside lever area
				if (e->x() >= margin_ && e->x() <= width()-margin_ && e->y() >= height()-margin_-lever_area_height_ && e->y() <= height()-margin_)
				{
					SignedInt pos = SignedInt(100*(e->x()-margin_)/float(gradient_area_width_));
					//be careful not to remove other levers...
					if (pos!=selected_ && !gradient_.exists(pos))
					{
						gradient_.remove(selected_);
						gradient_.insert(pos,selected_color_);
						selected_ = pos;
						repaint();
					}
				}
			}
		}
	
	void MultiGradientSelector::mouseReleaseEvent ( QMouseEvent * e )
	{
		if ( e->button() != Qt::LeftButton ) 
		{
			e->ignore();
			return;
	  } 	
		left_button_pressed_ = false;
	}
	
	void MultiGradientSelector::mouseDoubleClickEvent ( QMouseEvent * e )
	{
		for (UnsignedInt i=0;i<gradient_.size();++i)
		{
			SignedInt pos = SignedInt(float(gradient_.position(i))/100.0*gradient_area_width_+margin_+1);
			if (e->x() >= pos-3 && e->x() <= pos+4 && e->y() >= height()-margin_-lever_area_height_+8 && e->y() <= height()-margin_-lever_area_height_+15)
			{
				gradient_.insert(gradient_.position(i),QColorDialog::getColor(gradient_.color(i),this));
				if (SignedInt(gradient_.position(i))==selected_)
				{
					selected_color_=gradient_.color(i);
				}
				return;
			}
		}
	}
	
	void MultiGradientSelector::keyPressEvent ( QKeyEvent * e )
	{
		if (e->key() == Qt::Key_Delete)
		{
			gradient_.remove(selected_);
			selected_=-1;
			selected_color_=Qt::white;
			repaint();
		}
		else
		{
			e->ignore();		
		}
	}
	
	void MultiGradientSelector::stairsInterpolation(bool state)
	{
		if (state)
		{
			gradient_.setInterpolationMode(MultiGradient::IM_STAIRS);
		}
		else
		{
			gradient_.setInterpolationMode(MultiGradient::IM_LINEAR);
		}
	}
	
	void MultiGradientSelector::setInterpolationMode(MultiGradient::InterpolationMode mode)
	{
		gradient_.setInterpolationMode(mode);
	}
	
	UnsignedInt MultiGradientSelector::getInterpolationMode() const
	{
		return gradient_.getInterpolationMode();
	}

} //namespace OpenMS



