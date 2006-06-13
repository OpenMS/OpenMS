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
// $Id: WindowMower.C,v 1.5 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>

using namespace std;
namespace OpenMS
{

  //const String WindowMower::info_ = "augments the [peakcount] highest peaks in windows of size [windowsize]";
  
  WindowMower::WindowMower()
    : PreprocessingFunctor()
  {
		name_ = WindowMower::getName();
    defaults_.setValue("windowsize", 50); // smallest amino acid
    defaults_.setValue("peakcount", 2); // b and y ion
		param_ = defaults_;
  }

  WindowMower::WindowMower(const WindowMower& source)
    : PreprocessingFunctor(source)
  {
  }

  WindowMower& WindowMower::operator=(const WindowMower& source)
  {
    PreprocessingFunctor::operator=(source);
    return *this;
  }

  WindowMower::~WindowMower()
  {
  }
/*
  void WindowMower::operator()(MSSpectrum< DPeak<1> >& spec) const
  {
    double windowsize = (double)param_.getValue("windowsize");
    uint peakcount = (int)param_.getValue("peakcount");
    map<double,double> peaksinwindow; // peakheight,pos
    map<double,int> marks; // peaks get marked if they belong to the <peakcount> highest in the window
    for (MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end(); ++it)
    {
       peaksinwindow.clear();
       for (uint i = 0; (it+i) != spec.end() && (it+i)->getIntensity() < it->getIntensity()+windowsize ; ++i)
       {
         peaksinwindow.insert(make_pair((it+i)->getIntensity(),(it+i)->getPosition()[0]));
       }
       map<double,double>::reverse_iterator it2 = peaksinwindow.rbegin();
       for (uint i = 0; i < peakcount && i < peaksinwindow.size(); ++i)
       {
         marks[(it2++)->second]++;
       }
       //todo do something with the marking, maybe multiply the peaks with it
       for (uint i = 0; i < spec.size(); ++i)
       {
         spec.getContainer()[i].setIntensity( spec.getContainer()[i].getIntensity() + spec.getContainer()[i].getIntensity()*marks[spec.getContainer()[i].getPosition()[0]]);
       }
            
    }
  }
*/

	/*
  String WindowMower::info() const
  {
    return info_;
  }*/

}
