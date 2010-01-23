// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
#ifndef OPENMS_VISUAL_LISTEDITOR_H
#define OPENMS_VISUAL_LISTEDITOR_H


#include<OpenMS/DATASTRUCTURES/StringList.h>
#include <QtGui/QDialog>
#include<QtGui/QListWidget>
#include <QtGui/QItemDelegate>

	class QPushButton;
	
namespace OpenMS
{
	namespace Internal
	{
		class ListTable;
		class ListEditorDelegate;
	}
	
	/**
		@brief Editor for editing int, double and string lists (including output and input file lists)
	*/	
	class OPENMS_DLLAPI ListEditor
		:public QDialog
	{
		Q_OBJECT
		
		public:
			//types of lists
			enum Type
			{
				INT,
				FLOAT,
				STRING,
				OUTPUT_FILE,
				INPUT_FILE
			};	
		
			///Constructor	
			ListEditor(QWidget* parent = 0,QString title = "");
			///returns modified list
			StringList getList() const;
			///sets list (and its type)that will be modified by user
			void setList(const StringList& list, ListEditor::Type type);
			///set restrictions for list elements
			void setListRestrictions(const String& restrictions);
			///set name of type
			void setTypeName(QString name);

		private:
			///List type
			ListEditor::Type type_;
			///displays the list
			Internal::ListTable *listTable_;
			///Delegate between view and model
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

	/**
		@brief Namespace used to hide implementation details from users.
		
	*/	
	namespace Internal
	{
		class OPENMS_DLLAPI ListTable
			: public QListWidget
		{
			Q_OBJECT
	
			public:
	
				//Default Constructor
				ListTable(QWidget* parent =0);
	
				//returns a list_
				StringList getList();
	
				//sets new list
				void setList(const StringList& list, ListEditor::Type type);

			public slots:
				void createNewRow();
				void removeCurrentRow();

			private:
				///List type
				ListEditor::Type type_;
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
			  void updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex& index) const;
				
				//sets Type of List
				void setType(const ListEditor::Type type);
				//sets restrictions for listelements
				void setRestrictions(const String& restrictions);
				///set name of type
				void setTypeName(QString name);
				///sets the fileName
				void setFileName(QString name);

			private:
				/// Not implemented => private
				ListEditorDelegate();
				///List type
				ListEditor::Type type_;
				///restrictions for list elements
				String restrictions_;
				///type name. used to distinguish output/input from string lists
				QString typeName_;
				///used to set input and output values in setModelData
				mutable QString file_name_;

				
		};
	}





} // namespace OpenMS
#endif //OPENMS_VISUAL_LISTEDITOR_H
