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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_LISTSTACK_H
#define OPENMS_VISUAL_LISTSTACK_H

#include <OpenMS/CONCEPT/Types.h>

//QT
#include <qlistview.h>
#include <qwidgetstack.h>

//STL
#include <string>
#include <map>
#include <vector>

namespace OpenMS 
{
	/**
		@brief Tree view combined with a widget stack.
		
		Displays and manages a tree view of itemsand a stack of widgets.
		The shown stack item is determined by the activated entry in the tree view.
		
		\image html ListStack.png
		
		In the above example image a ListStack is shown that consists of the tree view (left red rectangle)
		and the widget stack (right red triangle).
		
		@ingroup Visual
	*/
	class ListStack : public QWidget
	{
		Q_OBJECT

		public:
			///Constructor
			ListStack( QWidget * parent = 0, const char * name = 0 );
			///Destructor
			~ListStack();

			///Expands all nodes (subnodes are inserted unexpanded by default).
			void expand();

			/**
				@brief Adds a widget with a certain name to the stack.
			
				Creator and parent are needed to locate the position where to insert the widget.
			*/
			void addWidget(std::string name, QWidget* widget, void* creator, void* parent=0);

			///returns a pointer to the active widget
			QWidget* activeWidget();


		protected:
			QWidgetStack* stack_;
			QListView* list_;
			QListViewItem* last_;
			std::map<void*,QListViewItem*> w_to_item_;
			
		protected slots:
			void raiseWidget_( QListViewItem* ptr );

	};

}
#endif // OPENMS_VISUAL_LISTSTACK_H

