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
#include <OpenMS/VISUAL/EnhancedTabBar.h>
#include <QtGui/QMouseEvent>

using namespace std;

namespace OpenMS
{

	EnhancedTabBar::EnhancedTabBar( QWidget * parent) 
		: QTabBar(parent)
	{
		
	}
	
	EnhancedTabBar::~EnhancedTabBar()
	{
		
	}


	void EnhancedTabBar::mouseDoubleClickEvent ( QMouseEvent * e )
	{
		if ( e->button() != Qt::LeftButton ) 
		{
			e->ignore();
			return;
    }
    for (int i=0; i<this->count(); ++i)
    {
			if (tabRect(i).contains(e->pos()))
			{
				emit doubleClicked(i);
				break;
			}
		}
	}
	
} //namespace OpenMS	

