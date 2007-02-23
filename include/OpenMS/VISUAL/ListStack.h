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


#ifndef OPENMS_VISUAL_LISTSTACK_H
#define OPENMS_VISUAL_LISTSTACK_H

//QT
#include <QtGui/QTreeWidget>
#include <QtGui/QStackedWidget>

namespace OpenMS 
{
	/**
		@brief QTreeWidget combined with a QStackedWidget.
		
		Displays and manages a tree view of items and a stack of widgets.
		
		\image html ListStack.png
		
		In the above example image a ListStack is shown that consists of the tree view (left red rectangle)
		and the widget stack (right red rectangle).
		
		@ingroup Visual
	*/
	class ListStack 
		: public QWidget
	{
		Q_OBJECT

		public:
			///Constructor
			ListStack( QWidget * parent = 0);
			///Destructor
			~ListStack();

			///Expands all nodes (subnodes are inserted unexpanded by default).
			void expand();

			/**
				@brief Adds an entry to the stack.
			
				@param name The name displayed in the the QTreeWidget
				@param widget The widget to associate with the name
				@param creator A pointer to the owner of the page
				@param parent A pointer to the owner's parent widget
				@param highlight Activates this page if true
			*/
			void addWidget(std::string name, QWidget* widget, void* creator, bool highlight, void* parent=0);

			///returns a pointer to the active widget
			QWidget* activeWidget();

		protected:
			/// Widget stack
			QStackedWidget* stack_;
			/// Tree view
			QTreeWidget* tree_;
			/// The last inserted item
			QTreeWidgetItem* last_;
			/// Connection map between owners and TreeWidgetItems (for finding the right place to add children)
			std::map<void*,QTreeWidgetItem*> w_to_item_;
			/// Connection map between TreeWidgetItem and index in the QStackedWidget
			std::map<QTreeWidgetItem*, int> item_to_index_;
			
		protected slots:
			void raiseActiveWidget_();

	};

}
#endif // OPENMS_VISUAL_LISTSTACK_H

