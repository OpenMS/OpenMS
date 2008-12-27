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

#include<OpenMS/VISUAL/ListEditor.h>
#include<OpenMS/DATASTRUCTURES/String.h>
//#include <OpenMS/DATASTRUCTURES/StringList.h>

//für DIALOG
#include<QtGui/QPushButton>
#include<QtGui/QVBoxLayout>
#include<QtGui/QHBoxLayout>
#include<QtGui/QHeaderView>
#include <QtGui/QMessageBox>
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>
#include <vector>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
		ListEditorDelegate::ListEditorDelegate(QObject* parent)
			:QItemDelegate(parent)
		{
		}
		
		QWidget *ListEditorDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem& , const QModelIndex& index) const
		{
			if(type_ == ListTable::STRING && restrictions_ != "")
			{
					QComboBox* editor = new QComboBox(parent);
					QStringList list;
					list.append("");
					list += restrictions_.toQString().split(",");
					editor->addItems(list);	
					return editor;
			}
			else
			{
				QLineEdit *editor = new QLineEdit(parent);
				editor->setFocusPolicy(Qt::StrongFocus);
				return editor;
			}	
		}
		
		void ListEditorDelegate::setEditorData(QWidget* editor, const QModelIndex& index) const
		{
			if(index.isValid())
			{
				QString str = index.data(Qt::DisplayRole).toString();
				if(qobject_cast<QLineEdit*>(editor)==0 )
				{
					int index = static_cast<QComboBox*>(editor)->findText(str);
					if (index==-1)
					{
						index = 0;
					}
					static_cast<QComboBox*>(editor)->setCurrentIndex(index);	
				}
				else
				{
					static_cast<QLineEdit*>(editor)->setText(str);
				}
			}
		}
		
		void ListEditorDelegate::setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const
		{
			QVariant present_value = index.data(Qt::DisplayRole);
			QVariant new_value;
			if(index.column()==0)
			{
				if(qobject_cast<QLineEdit*>(editor)==0 ) //Drop-down list for enums
				{
					new_value = QVariant(static_cast<QComboBox*>(editor)->currentText());
				}
				else
				{
					new_value = QVariant(static_cast<QLineEdit*>(editor)->text());
				}
				//check if it matches the restrictions or is empty
				if (new_value.toString()!="")
				{
					bool restrictions_met = true;
					switch(type_)
					{
						//check if valid integer
						case ListTable::INT:
						{
							bool ok;
							new_value.toString().toLong(&ok);
							if (!ok)
							{
								QMessageBox::warning(0,"Invalid value",QString("Cannot convert '%1' to integer number!").arg(new_value.toString()) );
								new_value = present_value;
							}
							//restrictions
							vector<String> parts;
							if (restrictions_.split(' ',parts))
							{
								if (parts[0]!="" && new_value.toInt()<parts[0].toInt())
								{
									restrictions_met = false;
								}
								if (parts[1]!="" && new_value.toInt()>parts[1].toInt())
								{
									restrictions_met = false;
								}
							}
						}
						break;
						case ListTable::FLOAT: //check if valid float
						{
							bool ok;
							new_value.toString().toDouble(&ok);
							if (!ok)
							{
								QMessageBox::warning(0,"Invalid value",QString("Cannot convert '%1' to floating point number!").arg(new_value.toString()) );
								new_value = present_value;
							}
							//restrictions
							vector<String> parts;
							if (restrictions_.split(' ',parts))
							{
								if (parts[0]!="" && new_value.toDouble()<parts[0].toDouble())
								{
									restrictions_met = false;
								}
								if (parts[1]!="" && new_value.toDouble()>parts[1].toDouble())
								{
									restrictions_met = false;
								}
							}
						}
						break;
						default:
						{
						}
					}	
					if(!restrictions_met)
					{
						QMessageBox::warning(0,"Invalid value",QString("Value restrictions not met: %1").arg(index.sibling(index.row(),3).data(Qt::DisplayRole).toString()) );
						new_value = present_value;
					}
				}
			}
			

				//check if modified
				if(new_value!=present_value)
				{
					model->setData(index, new_value);
					model->setData(index,QBrush(Qt::yellow),Qt::BackgroundRole);
					//emit modified(true);
				}
		}
		
		 
		void ListEditorDelegate::updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex& ) const
		{
			editor->setGeometry(option.rect);
		}

		void ListEditorDelegate::setType(const Internal::ListTable::Type type)
		{
			type_ = type;
		}
		void ListEditorDelegate::setRestrictions(const String& restrictions)
		{
			restrictions_ = restrictions;
		}
	
//////////////////////////////////////////////////////////////
//LISTTABLE
//////////////////////////////////////////////////////////////	
	ListTable::ListTable(QWidget* parent)
	:QTableWidget(parent)
	{	
		if(horizontalHeader() != NULL)
		{
			horizontalHeader() -> hide();
		}
		if(verticalHeader() != NULL)
		{
			verticalHeader() ->hide();
		}
	}

	ListTable::ListTable(int rows, int columns, QWidget* parent)
		:QTableWidget(rows,columns,parent)
	{
			if(horizontalHeader() != NULL)
		{
			horizontalHeader() -> hide();
		}
		if(verticalHeader() != NULL)
		{
			verticalHeader() ->hide();
		}
	}

	StringList ListTable::getList()
	{
		String stringit;
		list_.clear();
		for(Int i = 0; i < rowCount(); ++i)
		{
			stringit = item(i,0)->text();
			if(stringit != "")
			{
				stringit.trim();
			}
			list_.push_back(stringit);
		}
		return list_;
	}

	void ListTable::setList(const StringList &list)
	{
		QTableWidgetItem* item = NULL;
		for(UInt i = 0; i < list.size(); ++i)
		{
			item = new QTableWidgetItem(list[i].toQString());
			item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
			createNewRow();
			setItem(i,0,item);
		}
		list_ = list;
	}

	void ListTable::createNewRow()
	{
		insertRow(rowCount());
		setCurrentCell(rowCount()-1,0);
		QTableWidgetItem *item = new QTableWidgetItem;
		item->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
		setItem(rowCount()-1,0,item);
	}

	void ListTable::removeCurrentRow()
	{
		removeRow(currentRow());
	}

}
////////////////////////////////////////////////////////////
//ListEditor
////////////////////////////////////////////////////////////
ListEditor::ListEditor(QWidget *parent)
	: QDialog(parent)
{
	listTable_ = new Internal::ListTable(0,1,this);
	listTable_-> setRowHidden(-1,true);
	listDelegate_ = new Internal::ListEditorDelegate(listTable_);
	listTable_->setItemDelegate(listDelegate_);
	
	removeRowButton_ = new QPushButton(tr("&delete Row"));
	newRowButton_ = new QPushButton(tr("&new Row"));
	newRowButton_->setDefault(true);
	OkButton_ = new QPushButton(tr("&Ok"));
	CancelButton_ = new QPushButton(tr("&Cancel"));

	connect(OkButton_,SIGNAL(clicked()),
					this,SLOT(accept()));
	connect(CancelButton_,SIGNAL(clicked()),
					this, SLOT(reject()));
	
	connect(newRowButton_, SIGNAL(clicked()),
					listTable_, SLOT(createNewRow()));
	connect(removeRowButton_, SIGNAL(clicked()),
					listTable_, SLOT(removeCurrentRow()));
	
	QVBoxLayout *rightLayout = new QVBoxLayout;
	rightLayout -> addWidget(newRowButton_);
	rightLayout->addWidget(removeRowButton_);
	rightLayout->addWidget(OkButton_);
	rightLayout->addWidget(CancelButton_);
	rightLayout->addStretch();
	
	QHBoxLayout *mainLayout = new QHBoxLayout;
	mainLayout->addWidget(listTable_);
	mainLayout->addLayout(rightLayout);
	setLayout(mainLayout);
	
	setWindowTitle(tr("List Editor"));
}

StringList ListEditor::getList() const
{
	return listTable_->getList();
}

void ListEditor::setList(const StringList& list, Internal::ListTable::Type type)
{
	listTable_->setList(list);
	listDelegate_->setType(type);
}
	void ListEditor::setListRestrictions(const String& restrictions)
{
	listDelegate_->setRestrictions(restrictions);
}

}//namespace OpenMS
