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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_PARAMEDITOR_H
#define OPENMS_VISUAL_PARAMEDITOR_H

#include <OpenMS/CONCEPT/Types.h>


#include <QtGui/QItemDelegate>
#include <QtGui/QTreeWidget>

class QModelIndex;
class QStyleOptionViewItem;
class QAbstractItemModel;
class QStringList;
class QString;

namespace OpenMS
{
	class String;
	class Param;
	
	/**
		@brief Namespace used to hide implementation details from users.
		
	*/	
	namespace Internal
	{
		/**
			@brief Internal delegate class for ParamEditor
				
			This handles editing of items.
		*/
		class ParamEditorDelegate 
			: public QItemDelegate
		{
			Q_OBJECT

		 	public:
		 		///Constructor
			  ParamEditorDelegate(QObject *parent = 0);
				/// Returns the widget(combobox or QLineEdit) used to edit the item specified by index for editing. Prevents edit operations on nodes' values and types
			  QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;
				/// Sets the data to be displayed and edited by the editor for the item specified by index.
			  void setEditorData(QWidget *editor, const QModelIndex &index) const;
				/// Sets the data for the specified model and item index from that supplied by the editor. If data changed in a cell, that is if it is different from an initial value, then set its background color to yellow and emit the modified signal otherwise make it white
			  void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const;
				/// Updates the editor for the item specified by index according to the style option given.    
			  void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const;
			
			signals:
				/// signal for showing ParamEditor if the Model data changed
				void modified(bool) const;
		};
	}
	
	/**
		@brief A GUI for editing or viewing a Param object

		@improvment When loosing the focus, edit mode should be left (Marc)
		@improvment Prevent items/sections with the same name (Marc)
		
		@ingroup Visual
	*/
	class ParamEditor  
		: public QTreeWidget
	{
		Q_OBJECT
		
		public:
			/// Role of the entry
			enum
			{
				NODE,				///< Section
				NORMAL_ITEM,	///< Item that is always shown
				ADVANCED_ITEM	///< Item that is shown only in advanced mode
			};

			/// constructor
			ParamEditor(QWidget* parent=0);
			/// load method for const Param object
			void load(const Param& param);
			/// load method for editable Param object
			void loadEditable(Param& param);
			/// store edited data in Param object
			void store();
			/// is data changed since last save?
			bool isModified();
			/// Creates default shortcuts for copy, cut, paste, ...
			void createShortcuts();
			
		signals:
			/// item was edited
			void modified(bool);
		
		public slots:
			/// Switches between normal and advanced mode
			void toggleAdvancedMode(bool advanced);
			
		protected slots:
			/// deletes an item and its children
			void deleteItem();
			/// inserts an item
			void insertItem();
			/// inserts a node
			void insertNode();
			/// copy subtree
			void copySubTree();
			/// paste subtree
			void pasteSubTree();
			/// cut subtree
			void cutSubTree();
			/// Notifies the widget that the content was changed.
			/// Emits the modified(bool) signal if the state changed.
			void setModified(bool is_modified);
			/// Toggles between normal and advanced parameter mode of the selected item
			void toggleItemMode();
		protected:
			/// used to insert or delete elements by mouseclick events
			void contextMenuEvent(QContextMenuEvent* event);
			/// recursive helper method for method storeRecursive()
			void storeRecursive_(QTreeWidgetItem* child, String path, std::map<String,String>& section_descriptions);
			/// recursive helper method for slot deleteItem()
			void deleteItemRecursive_(QTreeWidgetItem* item);
			
			/// Param object for load(const Param&)
			Param* param_editable_;
			/// Param object for loadEditable(Param&)
			const Param* param_const_;         
			/// selected item
			QTreeWidgetItem* selected_item_;
			/// copied item 
			QTreeWidgetItem* copied_item_;
			/// Indicates that the data was modified since last store/load operation
			UInt modified_;
			/// Indicates if normal mode or advanced mode is activated
			bool advanced_mode_;
			
	};


} // namespace OpenMS

#endif // OPENMS_VISUAL_PARAMEDITOR_H
