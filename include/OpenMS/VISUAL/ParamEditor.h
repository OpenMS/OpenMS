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
//
// --------------------------------------------------------------------------
// $Maintainer: Stefan Rink $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_PARAMEDITOR_H
#define OPENMS_VISUAL_PARAMEDITOR_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QTreeWidget>
#include <string>

namespace OpenMS
{
 	class Param;
	
	/**
	
		@brief Visualization of class Param's loaded XML-file
		
		This class provides visualization for the XML-files and class Param.
		It can also be used to edit the data visually.
		
		@todo Allow valid Param types only (Stefan)
		@todo Allow deleting an entry (Stefan)
		@todo Show only entries below '1' (Stefan)
		@todo Make expanding/collapsing only a subtree possible (Stefan)
		@todo Add INIFileEditor TOPP tool (Stefan)
		
		@ingroup Visual
	*/
	class ParamEditor  
		: public QTreeWidget
	{
		Q_OBJECT
		
		public:
			/// constructor
			ParamEditor(QWidget* parent=0);
			/// load method for const Param object
			void load(const Param& param);
			/// load method for editable Param object
			void loadEditable(Param& param);
			/// used to insert or delete elements by certain key events
			void keyPressEvent(QKeyEvent* e);
			/// used to insert or delete elements by mouseclick events
			void contextMenuEvent(QContextMenuEvent* event);
			/// store edited data in Param object
			bool store() const;
			/// check if edited data still valid before storing
			bool isValid() const;
			/// delete all items
			void deleteAll();

		public slots:
			/// deletes an item and its children
			void deleteItem();
			/// inserts an item
			void insertItem();
			
		private:
			/// recursive helper method for method isValid()
			bool isValidRecursive_(QTreeWidgetItem* parent) const;
			/// recursive helper method for method storeRecursive()
			void storeRecursive_(const QTreeWidgetItem* child, String path) const;
			/// recursive helper method for slot deleteItem()
			void deleteItemRecursive_(QTreeWidgetItem* item);
			
		protected:
			/// Param object for load(const Param&)
			Param* param_editable_;
			/// Param object for load_editable(Param&)
			const Param* param_const_;         
			
	};


} // namespace OpenMS

#endif // OPENMS_VISUAL_PARAMEDITOR_H
