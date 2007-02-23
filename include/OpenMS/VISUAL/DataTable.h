// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: stefan_heess  $
// --------------------------------------------------------------------------

 
#ifndef OPENMS_VISUAL_DATATABLE_H
#define OPENMS_VISUAL_DATATABLE_H

#include <OpenMS/CONCEPT/Types.h>

//QT
#include <QtGui/QWidget>

class QPushButton;
class QGridLayout;
class QLineEdit;
class QTextEdit;
class QComboBox;

namespace OpenMS 
{
	/**
		@brief A class that provides some functions for displaying data.
		
		This class is a basic class for all classes to be displayed in the MetaData viewer. 
		So it provides some functions needed in all subclasses.
	*/
	class DataTable 
		: public QWidget
	{
		Q_OBJECT

		public: 
			/// Default constructor
			DataTable(bool editable, QWidget *parent =0);
			/// Adds a label to the grid layout.
			void addLabel(const QString &label);
		  /// Adds a line edit field with label to the grid layout.
			void addLineEdit(QLineEdit* &ptr ,  const QString &label);
			/// Adds a line edit field with label and button to the next free position in the grid.
			void addLineEditButton(const QString &labelname, QLineEdit* &ptr1, QPushButton* &ptr2, const QString &buttonlabel);
	
			/// Adds a text edit field to the grid layout.
			void addTextEdit(QTextEdit* &ptr ,  const QString &label);
			/// Adds a line edit field to the grid layout.
			void addComboBox(QComboBox* &ptr ,  const QString &label);//,  const std::string* items[]);
			/// Fills a combo box with data.
			void fillComboBox(QComboBox* &ptr , const std::string* items, int agr);
					
			/// Adds vertical spacer.
			void addVSpacer();
			/// Adds a button to the next free position in the grid.
			void addButton(QPushButton* &ptr, const QString &label);
			/// Adds two buttons in a row.
			void addHorizontalButtons(QPushButton* &ptr1, const QString &label1, QPushButton* &ptr2, const QString &label2);
			
			/// Adds a horizontal line as a seperator.
			void addSeperator();
			
			/// Adds an empty line.
			void addEmptyLine();
			
			/// Returns if the values are editable
			bool isEditable() const;
		
		
		protected:
		  /// The main layout.
			QGridLayout* mainlayout_;	
			
			/// Counter for the grid row.
			UnsignedInt row_;
			/// Counter for the grid column.
			UnsignedInt column_;
			
			/// Adds a label.
			void addLabel_(const QString &labelName, UnsignedInt row );

		private:
			/// Edit flag
			bool editable_;
				
	};


}
#endif
