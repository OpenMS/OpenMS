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
// $Id: FactoryProductView.C,v 1.2 2006/03/28 08:03:39 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
#include <OpenMS/VISUAL/FactoryProductView.h>

#include <sstream>
#include <map>

using namespace std;

namespace OpenMS
{

FactoryProductView::FactoryProductView(QWidget* parent, const char* name)
  :QLabel(parent,name)
{
}

FactoryProductView::~FactoryProductView()
{
}

void FactoryProductView::displayFactoryProduct(const FactoryProduct* const conf)
{
  ostringstream ss;
  ss << "<qt>";
  ss << "<u>Name:</u> \t" << conf->getName() << "<br>"
   /* << "<u>Info:</u> \t" << conf->info() << "<br>"*/;
  if (conf->getParam().size())
  {
    for (Param::ConstIterator cmit = conf->getParam().begin(); cmit != conf->getParam().end(); ++cmit)
    {
     ss << "<u>" << cmit->first << ":</u>\t" << (string)cmit->second << "<br>";
    }
  }
  else
  {
    ss << "none\n";
  }
  ss << "</qt>";
  setText(ss.str().c_str());
}

}
