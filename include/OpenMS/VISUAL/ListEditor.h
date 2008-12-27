// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: David Wojnar $
// --------------------------------------------------------------------------
#ifndef LIST_EDITOR_H
#define LIST_EDITOR_H


#include<OpenMS/DATASTRUCTURES/StringList.h>
#include <QtGui/QDialog>
#include<QtGui/QTableWidget>
#include <QtGui/QItemDelegate>

	class QPushButton;
	
namespace OpenMS
{

	/**
		@brief Namespace used to hide implementation details from users.
		
	*/		
		namespace Internal
	{
			class ListTable
			: public QTableWidget
			{
			Q_OBJECT
	
			public:
			//types of lists
				enum Type
				{
					EMPTY_VALUE,
					INT,
					FLOAT,
					STRING,
					OUTPUT_FILE,
					INPUT_FILE
				};
	
				//Default Constructor
				ListTable(QWidget* parent =0);
	
				ListTable(int rows, int columns, QWidget* parent = 0);
	
				//returns a list_
				StringList getList();
	
				//sets new list
				void setList(const StringList& list);

			public slots:
				void createNewRow();
				void removeCurrentRow();

			private:
			//everything is internally stored as stringlist
				StringList list_;
		};	
		
		/**
			@brief Internal delegate class
			
			This handles editing of items.
		*/
		class OPENMS_DLLAPI ListEditorDelegate 
			: public QItemDelegate
		{
			Q_OBJECT

		 	public:
		 		///Constructor
			  ListEditorDelegate(QObject* parent);
				/// not reimplemented
			  QWidget *createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const;
				/// Sets the data to be displayed and edited by the editor for the item specified by index.
			  void setEditorData(QWidget* editor, const QModelIndex& index) const;
				/// Sets the data for the specified model and item index from that supplied by the editor. If data changed in a cell, that is if it is different from an initial value, then set its background color to yellow and emit the modified signal otherwise make it white
			  void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const;
				/// Updates the editor for the item specified by index according to the style option given.    
			  void updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex &index) const;
				
				//sets Type of List
				void setType(const ListTable::Type type);
				//sets restrictions for listelements
				void setRestrictions(const String& restrictions);

			private:
				/// Not implemented
				ListEditorDelegate();
				///List type
				Internal::ListTable::Type type_;
				///restrictions for list elements
				String restrictions_;
		};
	}
//DIALOG
class ListEditor
:public QDialog
{
	Q_OBJECT
	
	public:
		///Constructor	
		ListEditor(QWidget* parent = 0);
		///returns modified list
		StringList getList() const;
		///sets list (and its type)that will be modified by user
		void setList(const StringList& list, Internal::ListTable::Type type);
		///set restrictions for list elements
		void setListRestrictions(const String& restrictions);
		/// get type of list
		Internal::ListTable::Type getType();
	
	private:
		///Inherits QTableWidget
		Internal::ListTable *listTable_;
		///Inherits QItemDelegate
		Internal::ListEditorDelegate *listDelegate_;
		/// buttton for new Row
		QPushButton *newRowButton_;
		///button for removing row
		QPushButton *removeRowButton_;
		///button clicked if modifications are accepted
		QPushButton *OkButton_;
		///button clicked if modifications are rejected
		QPushButton *CancelButton_;
};




} // namespace OpenMS
#endif //LIST_EDITOR_H
