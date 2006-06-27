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
// $Id: MultiGradient.C,v 1.8 2006/03/28 08:03:39 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/MultiGradient.h>

// stl
#include <iostream>
#include <sstream>
#include <math.h>

using namespace std;
using namespace OpenMS;
using namespace OpenMS::Exception;

	MultiGradient::MultiGradient():pos_col_(), interpolation_mode_(IM_LINEAR)
	{
		pos_col_[0] = Qt::white;
		pos_col_[100] = Qt::black;		
	}
	
	MultiGradient::~MultiGradient()
	{
		
	}

	UnsignedInt MultiGradient::size() const
	{
		return pos_col_.size();
	}
	

	void MultiGradient::insert (SignedInt position, const QColor& color)
	{
		if (position >= 0 || position <=100 )
		{
			pos_col_[position]=color;
		}
	}

	bool MultiGradient::remove (SignedInt position)
	{
		if (position < 1 || position > 99 )
		{
			return false;
		}
		
		map<UnsignedInt,QColor>::iterator it = pos_col_.find(position);
		if (it != pos_col_.end())
		{
			pos_col_.erase(it);
			return true; 
		}
		return false;
	}


	UnsignedInt MultiGradient::position(UnsignedInt index) throw(IndexUnderflow,IndexOverflow)
	{
		if (index>size()-1)
		{
			throw IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}

		map<UnsignedInt,QColor>::iterator it = pos_col_.begin();
		for (UnsignedInt i=0; i<index; ++i)
		{
			++it;
		}		
		return it->first;
	}
	
	const QColor& MultiGradient::color(UnsignedInt index) throw(IndexUnderflow,IndexOverflow)
	{
		if (index>size()-1)
		{
			throw IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		
		map<UnsignedInt,QColor>::iterator it = pos_col_.begin();
		for (UnsignedInt i=0; i<index; ++i)
		{
			++it;
		}
		return it->second;
	}

	QColor MultiGradient::interpolatedColorAt(double position) const
	{
		if  (position <= 0.0 )
		{
			return pos_col_.begin()->second;
		}
		
		if  (position >= 100.0 )
		{
			return (--(pos_col_.end()))->second;
		}
		
		//linear
		if (interpolation_mode_==IM_LINEAR)
		{
			map<UnsignedInt,QColor>::const_iterator it1 = pos_col_.lower_bound(SignedInt(position));
			if (it1->first == UnsignedInt(position))
			{
				return it1->second;
			}
			else
			{
				map<UnsignedInt,QColor>::const_iterator it0 = it1;
				--it0;
				double factor = (position-it0->first)/(it1->first-it0->first);
				return QColor(SignedInt(factor*it1->second.red()+(1-factor)*it0->second.red()) 
				            , SignedInt(factor*it1->second.green()+(1-factor)*it0->second.green()) 
				            , SignedInt(factor*it1->second.blue()+(1-factor)*it0->second.blue())  );
			}
		}
		//stairs
		else
		{
			map<UnsignedInt,QColor>::const_iterator it = pos_col_.upper_bound(SignedInt(position));
			--it;
			return it->second;
		}
		
		
	}

	QColor MultiGradient::interpolatedColorAt(double position, double min, double max) const
	{
		return interpolatedColorAt((position-min)/(max-min)*100);
	}

	void MultiGradient::setInterpolationMode(UnsignedInt mode)
	{
		if (mode == IM_LINEAR || mode == IM_STAIRS)
		{
			interpolation_mode_ = mode;
		}
	}

	UnsignedInt MultiGradient::getInterpolationMode() const
	{
		return interpolation_mode_;
	}


	string MultiGradient::toString() const
	{
		stringstream out;
		
		//Interpolation Mode
		if (getInterpolationMode()==IM_LINEAR)
		{
			out << "Linear|";
		}
		else if (getInterpolationMode()==IM_STAIRS)
		{
			out << "Stairs|";
		}		
		
		for (map<UnsignedInt,QColor>::const_iterator it = pos_col_.begin(); it!=pos_col_.end(); ++it )
		{
			if (it!=pos_col_.begin())
			{
				out << ";";
			}
			out << it->first << "," << it->second.name();
		}
		return out.str();
	}
	
	void MultiGradient::fromString(const string& gradient)
	{
		pos_col_.clear();

		if (gradient == "")
		{
			pos_col_[0] = Qt::white;
			pos_col_[100] = Qt::black;			
			return;		
		}

		string g(gradient);
		string::iterator tmp(g.begin());
		UnsignedInt tmp_pos;
		for (string::iterator it = g.begin(); it!=(g.end()+1);++it)
		{
			if (*it == '|')
			{
				//interploation mode
				if (string(tmp,it)=="Linear")
				{
					setInterpolationMode(IM_LINEAR);
				}
				else if (string(tmp,it)=="Stairs")
				{
					setInterpolationMode(IM_STAIRS);
				}
			}
			else if (*it == ';' || it == g.end())
			{
				pos_col_[tmp_pos] = QColor(string(tmp,it).c_str());
				tmp = it+1;
			}
			else if (*it == ',')
			{
				tmp_pos = atoi(string(tmp,it).c_str());
				tmp = it+1;				
			}
		}
	}
	
	void MultiGradient::activatePrecalculationMode(double min, double max, UnsignedInt steps)
	{
		pre_.clear();
		double step_length = fabs(max - min)/steps;
		pre_[min] = interpolatedColorAt(0);
		for (UnsignedInt step = 1; step < steps; ++step)
		{
			//cout << step << " -> "<<min+step_length*step<<endl;
			pre_[min+step_length*step] = interpolatedColorAt(double(step)/steps*100);
		}
		pre_[max] = interpolatedColorAt(100);
	}
	
	void MultiGradient::deactivatePrecalculationMode()
	{
		pre_.clear();
	}

	const QColor& MultiGradient::precalculatedColorAt(double position) const throw (Exception::OutOfSpecifiedRange)
	{
		if (pre_.size()==0)
		{
			throw OutOfSpecifiedRange(__FILE__, __LINE__, __PRETTY_FUNCTION__, position,0,0);
		}
			
		//lookup
		map<double,QColor>::const_iterator it = pre_.lower_bound(position);
		// error handling
		if (it==pre_.end())
		{
			--it;
		}
		//return ref
		return it->second;
	}

	bool MultiGradient::exists (SignedInt position)
	{
		return pos_col_.find(position)!=pos_col_.end();
	}

