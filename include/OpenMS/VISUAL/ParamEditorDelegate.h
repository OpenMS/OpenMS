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

#ifndef OPENMS_VISUAL_PARAMEDITORDELEGATE_H
#define OPENMS_VISUAL_PARAMEDITORDELEGATE_H

#include <OpenMS/CONCEPT/Types.h>
#include <QtGui/QItemDelegate>

class QModelIndex;
class QStyleOptionViewItem;
class QAbstractItemModel;
class QStringList;
namespace OpenMS
{
	/**
		
		@brief Delegate class for ParamEditor
			
		This class provides modifies visualization for the items in ParamEditor.
		It places a Combobox in the second column and prevents edit operations on nodes' values and types
			
			
		@ingroup Visual
	*/
	class ParamEditorDelegate : public QItemDelegate
	 {
	     Q_OBJECT

	 public:
	     ParamEditorDelegate(QObject *parent = 0);

	     QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
				   const QModelIndex &index) const;

	     void setEditorData(QWidget *editor, const QModelIndex &index) const;
	     void setModelData(QWidget *editor, QAbstractItemModel *model,
			       const QModelIndex &index) const;

	     void updateEditorGeometry(QWidget *editor,
		 const QStyleOptionViewItem &option, const QModelIndex &index) const;
	 private:
		 QStringList combo_list_;
	 };
 }
#endif

