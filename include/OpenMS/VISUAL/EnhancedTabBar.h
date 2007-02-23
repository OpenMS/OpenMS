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

#ifndef OPENMS_VISUAL_ENHANCEDTABBAR_H
#define OPENMS_VISUAL_ENHANCEDTABBAR_H

#include <OpenMS/CONCEPT/Types.h>

//QT
#include <QtGui/QTabBar>
class QMouseEvent;

namespace OpenMS 
{
	/**
		@brief Tab bar which is aware of double clicking.
		
		It emits a signal, when a tab is double clicked.
		
		@ingroup Visual
	*/
	class EnhancedTabBar
		: public QTabBar
	{
		Q_OBJECT
		public:
		/// Constructor
		EnhancedTabBar( QWidget * parent = 0);
		/// Destructor
		~EnhancedTabBar();
	
		signals:
		/// Signal emited when double clicked. Returns the tab index
		void doubleClicked(int);
	
		protected:
		void mouseDoubleClickEvent ( QMouseEvent * e );
	};

}
#endif // OPENMS_VISUAL_ENHANCEDTABBAR_H

