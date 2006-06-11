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
// $Id: ResultView.C,v 1.3 2006/03/28 08:03:39 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
#include <OpenMS/VISUAL/ResultView.h>

#include <sstream>
#include <iostream>

using namespace std;

namespace OpenMS
{
    
  ResultView::ResultView(QWidget* parent , const char* name )
    :QLabel(parent,name)
  {
  }

  ResultView::~ResultView()
  {
  }

  void ResultView::displayResults(const map<string,double>& result)
  {
    ostringstream ss;
    ss << "<qt> <font color=\"#ff0000\">";
    for (map<string,double>::const_iterator cmit = result.begin(); cmit != result.end(); ++cmit)
    {
      ss << cmit->first << "\t" << cmit->second << "<br>";
    }
    ss << "</font></qt>";
    setText(ss.str().c_str());
  }

}
