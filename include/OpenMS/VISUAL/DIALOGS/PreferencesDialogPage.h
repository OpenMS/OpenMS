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


#ifndef OPENMS_VISUAL_DIALOGS_PREFERENCESDIALOGPAGE_H
#define OPENMS_VISUAL_DIALOGS_PREFERENCESDIALOGPAGE_H

// QT
#include <QtGui/QWidget>
class QGridLayout;
class QString;
class QGroupBox;
class QLayout;
class QSpinBox;

// STL
#include <string>

namespace OpenMS
{
	class PreferencesManager;
	
	/**
		@brief Base class for all PreferencesManager dialog pages.
		
		
		
		@ingroup Dialogs
	*/
	class PreferencesDialogPage
		: public QWidget
	{
		Q_OBJECT

		public:
			///constructor
			PreferencesDialogPage( PreferencesManager* manager, QWidget* parent = 0);
			///destructor
			virtual ~PreferencesDialogPage();

			/// returns the help text
			const std::string& getHelpText() const;
			///sets the help text
			void setHelpText(const std::string& text);

			/// load values from the PreferencesManager's Param
			virtual void load();
			/// write changes in Dialog to the PreferencesManager's Param
			virtual void save();



		protected:
			///Adds a QSpinBox to the @p parent
			QSpinBox* addSpinBox(QWidget* parent, int min, int max, int step) const;
			///Adds a new QGroupBox with a @p label to the QGridlayout @p grid
			QGroupBox* addBox(QGridLayout* grid, int row, int col, const QString& label, int row_span=1, int col_span=1);
			///Adds a @p label and a @p widget to the give @p row of a QGridlayout @p grid
			void addWidget(QLayout* grid, int row, const QString& label, QWidget* widget) const;
			///Adds a @p label and a @p layout to the give @p row of a QGridlayout @p grid
			void addLayout(QLayout* grid, int row, const QString& label, QLayout* layout) const;
			///Appends a row and column with strech @p stretch to a QGridlayout @p grid
			void finish(QLayout* grid, int stretch = 2) const;
			///Pointer to the PreferencesManager
			PreferencesManager* manager_;
			///stores the help text
			std::string help_;

	};
} // namespace OpenMS


#endif

