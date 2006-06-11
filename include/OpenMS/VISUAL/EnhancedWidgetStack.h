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
// $Id: EnhancedWidgetStack.h,v 1.5 2006/03/03 17:55:35 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_ENHANCEDWIDGETSTACK_H
#define OPENMS_VISUAL_ENHANCEDWIDGETSTACK_H

#include <OpenMS/CONCEPT/Types.h>

//QT
#include <qwidgetstack.h>
#include <qlistview.h>

namespace OpenMS
{
	/**
		@brief A QWidgetStack that accepts QListItem* as id for the stack item to raise.
	
		@ingroup Visual
	*/
	class EnhancedWidgetStack: public QWidgetStack
	{
		Q_OBJECT
		public:
			EnhancedWidgetStack( QWidget * parent = 0, const char * name = 0 );
	
			~EnhancedWidgetStack();
	
			PointerSizeInt addWidget ( QWidget * w, QListViewItem* ptr);

		public slots:
			void raiseWidget ( QListViewItem* ptr );
	};
}
#endif // OPENMS_VISUAL_ENHANCEDWIDGETSTACK_H

