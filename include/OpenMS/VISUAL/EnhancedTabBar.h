// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_ENHANCEDTABBAR_H
#define OPENMS_VISUAL_ENHANCEDTABBAR_H

#include <OpenMS/config.h>

//QT
#include <QtGui/QTabBar>
class QMouseEvent;
class QMimeData;

namespace OpenMS 
{
	class String;

	/**
		@brief Convenience tab bar implementation
		
		This tab bar differs in several ways form the QTabBar:
		- you can close a tab by double-clicking it or through its context menu.
		- it works based on tab identifiers (a fixed id stored in tab data) rather than on tab indices, which might
		  change when inserting or removing a tab.
		- it accepts all drag-and-drop actions and emits signals to handle them.
		
		@ingroup Visual
	*/
	class OPENMS_DLLAPI EnhancedTabBar
		: public QTabBar
	{
		Q_OBJECT
		public:
			/// Constructor
			EnhancedTabBar( QWidget * parent = 0);
			/// Destructor
			~EnhancedTabBar();
			
			/// Adds a new tab with the name @p text and the identifier @p id
			int addTab(const String& text, int id);
			/// Selects the tab with identifier @p id
			void setCurrentId(int id);
			
		public slots:
			/// Remove the tab with identifier @p id
			void removeId(int id);
			
		signals:
			/// Signal that indicates that the current tab changed
			void currentIdChanged(int id);
			/// Signal that indicates that the tab with identifier @p id is about to be removed (double click or context menu)
			void aboutToCloseId(int id);

			/// Signal that is emitted, when a drag-and-drop action ends on a tab
			void dropOnTab(const QMimeData* data, QWidget* source, int id);			
			/// Signal that is emitted, when a drag-and-drop action ends on the unused space on the right side of the tabs.
			void dropOnWidget(const QMimeData* data, QWidget* source);
			
		protected:
			///@name Remplemented Qt events
			//@{
			void mouseDoubleClickEvent(QMouseEvent* e);
			void contextMenuEvent(QContextMenuEvent* e);
			void dragEnterEvent(QDragEnterEvent* e);
			void dropEvent(QDropEvent* e);
			//@}
			
			///Returns the QTabBar index of the tab at position @p pos. If there is no tab at that position -1 is returned.
			int tabAt_(const QPoint& pos);

		protected slots:
			/// Slot that translates the currentChanged(int) signal to currentIdChanged(int)
			void currentChanged_(int id);
	};

}
#endif // OPENMS_VISUAL_ENHANCEDTABBAR_H

