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

#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <QtGui/QMessageBox>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>
#include <QtGui/QContextMenuEvent>
#include <QtGui/QShortcut>
#include <QtGui/QMenu>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
		ParamEditorDelegate::ParamEditorDelegate(QObject *parent)
			: QItemDelegate(parent)
		{
		}
		
		QWidget *ParamEditorDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem& /*option*/, const QModelIndex& index ) const
		{
			Int type = index.data(Qt::UserRole).toInt();
			// name -> QLineEdit with a regex validator that allows only words
			if(index.column()==0)	
			{
				QRegExp rx("[A-Za-z0-9_]*");
				QValidator *validator = new QRegExpValidator(rx, parent);
				QLineEdit *editor = new QLineEdit(parent);
				editor->setValidator(validator);
				editor->setFocusPolicy(Qt::StrongFocus);
				return editor;
			}
			// value -> QLineEdit
			else if(index.column()==1 && type==ParamEditor::ITEM)
			{
				QLineEdit *editor = new QLineEdit(parent);
				editor->setFocusPolicy(Qt::StrongFocus);
				return editor;
			}
			// type -> combobox to choose the type 
			else if (index.column() == 2 && type==ParamEditor::ITEM) 
			{
				QComboBox *editor = new QComboBox(parent);
				QStringList list;
				list<<"int"<<"float"<<"string";
				editor->addItems(list);
				editor->setFocusPolicy(Qt::StrongFocus);
				return editor;
			}
			// description -> QLineEdit
			else if(index.column()==3)
			{
				QLineEdit *editor = new QLineEdit(parent);
				editor->setFocusPolicy(Qt::StrongFocus);
				return editor;
			}
			return 0;
		}
		
		void ParamEditorDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
		{
			QString str = index.data(Qt::DisplayRole).toString();
			//name
			if(index.column()==0)
			{
				static_cast<QLineEdit*>(editor)->setText(str);
		 	}
		 	//value
			else if(index.column()==1)
			{
				static_cast<QLineEdit*>(editor)->setText(str);
			}
			//type
			else if(index.column()==2)
			{
				QStringList list;
				list<<"int"<<"float"<<"string";
				int pos = list.indexOf(str);
				if (pos!=-1)
				{
					static_cast<QComboBox*>(editor)->setCurrentIndex(pos);
				}
			}
			//description
			else if(index.column()==3)
			{
				static_cast<QLineEdit*>(editor)->setText(str);
			}
		}
		
		void ParamEditorDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
		{
			QVariant present_value = index.data(Qt::DisplayRole);
			QVariant new_value;
			//name
			if(index.column()==0)
			{
				new_value = QVariant(static_cast<QLineEdit*>(editor)->text());
				if (new_value.toString().isEmpty())
				{
					QMessageBox::warning(0,"Invalid name","A section name cannot be empty!");
					new_value = present_value;
				}
			}
			//value
			else if(index.column()==1)
			{
				QString type = index.sibling(index.row(),2).data(Qt::DisplayRole).toString();
				new_value = QVariant(static_cast<QLineEdit*>(editor)->text());
				if (type=="int") //check if valid integer
				{
					bool ok;
					new_value.toString().toLong(&ok);
					if (!ok)
					{
						QMessageBox::warning(0,"Invalid value",QString("Cannot convert '%1' to integer number!").arg(new_value.toString()) );
						new_value = present_value;
					}
				}
				else if (type=="float") //check if valid float
				{
					bool ok;
					new_value.toString().toDouble(&ok);
					if (!ok)
					{
						QMessageBox::warning(0,"Invalid value",QString("Cannot convert '%1' to floating point number!").arg(new_value.toString()) );
						new_value = present_value;
					}
				}
				
			}
			//type
			else if(index.column()==2)
			{
				new_value = QVariant(static_cast<QComboBox*>(editor)->currentText());
			}
			//description
			else if(index.column()==3)
			{
				new_value = QVariant(static_cast<QLineEdit*>(editor)->text());
			}
			
			//check if modified
			if(new_value!=present_value)
			{
				model->setData(index, new_value);
				model->setData(index,QBrush(Qt::yellow),Qt::BackgroundRole);
				emit modified(true);
			}
		}
		 
		void ParamEditorDelegate::updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex& /*index*/) const
		{
			editor->setGeometry(option.rect);
		}
	}

	ParamEditor::ParamEditor(QWidget * parent)
	  : QTreeWidget(parent),
	  	param_editable_(0),
	  	param_const_(0),
			selected_item_(0),
			copied_item_(0),
			modified_(false)
	{
		setMinimumSize(800,500);
		setItemDelegate(new Internal::ParamEditorDelegate);	// the delegate from above is set
		setWindowTitle("ParamEditor");
		setColumnCount(4);
		connect(itemDelegate(),SIGNAL(modified(bool)),this,SLOT(setModified(bool)));
		QStringList list;
		list.push_back("name");
		list.push_back("value");
		list.push_back("type");
		list.push_back("description");
		setHeaderLabels(list);
	}
	
	void ParamEditor::createShortcuts()
	{
		//cout << "creating shortcuts" << endl;
		new QShortcut(Qt::CTRL+Qt::Key_C, this, SLOT(copySubTree()));
		new QShortcut(Qt::CTRL+Qt::Key_X, this, SLOT(cutSubTree()));
		new QShortcut(Qt::CTRL+Qt::Key_V, this, SLOT(pasteSubTree()));
		new QShortcut(Qt::CTRL+Qt::Key_E, this, SLOT(insertNode()));
		new QShortcut(Qt::CTRL+Qt::Key_N, this, SLOT(insertItem())); 
		new QShortcut(Qt::Key_Delete, this, SLOT(deleteItem()));	
	}

	void ParamEditor::load(const Param& param)
	{
		clear();
		
		QTreeWidgetItem* parent=this->invisibleRootItem();
		param_const_=&param;
		QTreeWidgetItem* item = NULL;	
		
		for(Param::ParamIterator it=param.begin();it!=param.end();++it)
		{
			//********handle opened/closed nodes********
			const std::vector< Param::ParamIterator::TraceInfo >& trace = it.getTrace();
			for(std::vector< Param::ParamIterator::TraceInfo >::const_iterator it2 = trace.begin(); it2!=trace.end(); ++it2)
			{
				if (it2->opened) //opened node
				{
					item = new QTreeWidgetItem(parent);
					//name
					item->setText(0, it2->name.toQString());
					//description
					item->setText(3, it2->description.toQString());
					//role
					item->setData(0,Qt::UserRole,NODE);
					item->setData(1,Qt::UserRole,NODE);
					item->setData(2,Qt::UserRole,NODE);
					item->setData(3,Qt::UserRole,NODE);
					//flags
					if(param_editable_!=NULL)
					{
						item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
					}
					else 
					{
						item->setFlags( Qt::ItemIsEnabled );
					}
					parent=item;
				}
				else //closed node
				{
					parent=parent->parent();
					if(parent==NULL) parent=invisibleRootItem();
				}
			}
			
			//********handle item********
			item = new QTreeWidgetItem(parent);
			//handle user parameter
			if (it->user == true)
			{
				QFont font = item->font(0);
				font.setWeight(QFont::Bold);
				item->setFont(0, font);
				item->setFont(1, font);
				item->setFont(2, font);
				item->setFont(3, font);
			}
			//name
			item->setText(0, it->name.toQString());
			//value
			item->setText(1, QString::fromStdString(it->value.toString().c_str()));
			//type
			switch(it->value.valueType())
			{
				case DataValue::INTVALUE:
				case DataValue::LONVALUE:
				case DataValue::SHOVALUE:
					item->setText(2, "int");
					break;
				case DataValue::FLOVALUE:
				case DataValue::DOUVALUE:
					item->setText(2, "float");
					break;
				case DataValue::STRVALUE:
					item->setText(2, "string");
					break;
				default:
					break;
			};
			//description
			item->setText(3, it->description.toQString());
			//role
			item->setData(0,Qt::UserRole,ITEM);
			item->setData(1,Qt::UserRole,ITEM);
			item->setData(2,Qt::UserRole,ITEM);
			item->setData(3,Qt::UserRole,ITEM);
			//flags
			if(param_editable_!=NULL)
			{
				item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
			}
			else 
			{
				item->setFlags( Qt::ItemIsEnabled );
			}
		}

		expandAll();
		
		resizeColumnToContents(0);
		resizeColumnToContents(1);
		resizeColumnToContents(2);
		resizeColumnToContents(3);		
	}

	void ParamEditor::loadEditable(Param& param)
	{
		param_editable_=&param;
		load(param);
	}
	    
	void ParamEditor::store()
	{
		if(param_editable_!=NULL)
		{
			QTreeWidgetItem* parent=this->invisibleRootItem();
			param_editable_->clear();
		
			for (Int i = 0; i < parent->childCount();++i)
			{
				map<String,String> section_descriptions;
				storeRecursive_(parent->child(i),"", section_descriptions);	//whole tree recursively
			}	
		}
			
		setModified(false);
	}
	    
	void ParamEditor::deleteItem()
	{
		QTreeWidgetItem* item=selected_item_;
		selected_item_=NULL;
		if (!item)
		{
			item=currentItem();
			if(!item)
			{
				return;
			}
		}
		for(int i=item->childCount()-1; i>=0; i--)
		{
			deleteItemRecursive_(item->child(i));
		}
		delete item;
		setModified(true);
	}
		
	void ParamEditor::deleteItemRecursive_(QTreeWidgetItem* item)
	{				
		for(Int i=item->childCount()-1; i>=0; i--)	// childCount() is used every time and we count backwards because the number of childs changes
		{
			deleteItemRecursive_(item->child(i));
		}
		
		delete item;
	}
	
	void ParamEditor::deleteAll()
	{
		QTreeWidgetItem* item=invisibleRootItem();

		for (Int i = item->childCount()-1; i >=0;i--)
		{
			deleteItemRecursive_(item->child(i));
		}
		setModified(true);
	}
		
	void ParamEditor::insertItem()
	{
		if(!selected_item_)	// if no item is selected
		{
			selected_item_=currentItem();	// get the current item
			
			if(!selected_item_)	// if we didn't select any item get the root item
			{
				selected_item_=invisibleRootItem();
			}
		}
		
		QTreeWidgetItem* parent=selected_item_;
		if(parent->data(0,Qt::UserRole)==ITEM)	// if selected item is an item then get the parent node as parent for new item 
		{
			parent = selected_item_->parent();
		}
		QTreeWidgetItem* item = new QTreeWidgetItem(parent);
		item->setText(0, "name");
		item->setText(2, "string");
		item->setData(0,Qt::UserRole,ITEM);
		item->setData(1,Qt::UserRole,ITEM);
		item->setData(2,Qt::UserRole,ITEM);
		item->setData(3,Qt::UserRole,ITEM);
		item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
		item->setBackground (0,Qt::yellow);
		item->setBackground (1,Qt::yellow);
		item->setBackground (2,Qt::yellow);
		item->setBackground (3,Qt::yellow);
		setCurrentItem(item);
		editItem(item);
		setModified(true);
	}
	    
	void ParamEditor::insertNode()
	{
		if(!selected_item_)
		{
			selected_item_=currentItem();
			if(!selected_item_)
			{
				selected_item_=invisibleRootItem();
			}
		}
		
		QTreeWidgetItem* parent=selected_item_;
		if(parent->data(0,Qt::UserRole)==ITEM)
		{
			parent = selected_item_->parent();
		}
		QTreeWidgetItem* item = new QTreeWidgetItem(parent);
		item->setText(0, "name");
		item->setData(0,Qt::UserRole,NODE);
		item->setData(1,Qt::UserRole,NODE);
		item->setData(2,Qt::UserRole,NODE);
		item->setData(3,Qt::UserRole,NODE);
		item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
		item->setBackground (0,Qt::yellow);
		item->setBackground (1,Qt::yellow);
		item->setBackground (2,Qt::yellow);
		item->setBackground (3,Qt::yellow);
		setCurrentItem(item);
		editItem(item);
		setModified(true);
	}
	
	void ParamEditor::storeRecursive_(QTreeWidgetItem* child, String path,map<String,String>& section_descriptions)
	{
		child->setData ( 0, Qt::BackgroundRole, QBrush(Qt::white));
		child->setData ( 1, Qt::BackgroundRole, QBrush(Qt::white));
		child->setData ( 2, Qt::BackgroundRole, QBrush(Qt::white));
		child->setData ( 3, Qt::BackgroundRole, QBrush(Qt::white));
		
		if (path=="")
		{
			path = child->text(0).toStdString();
		}
		else
		{
			path = path + ":" + child->text(0).toStdString();	
		}
		
		String description = child->text(3).toStdString();
		
		if(child->text(2)=="") // node
		{
			if (description != "")
			{
				section_descriptions.insert(make_pair(path,description));
			}
		}
		else //item + section descriptions
		{			
			bool user = (child->font(0).weight() == QFont::Bold);
			if(child->text(2)=="float")
			{
				param_editable_->setValue(path,child->text(1).toDouble(),description,user);
			}
			else if(child->text(2)=="string")
			{
				param_editable_->setValue(path, child->text(1).toStdString(),description,user);
				//std::cerr<<"\n"<<path<<":  "<<child->text(1).toStdString()<<"\n";
			}
			else if(child->text(2)=="int")
			{
				param_editable_->setValue(path, child->text(1).toInt(),description,user);
			}

			// set description node description if the prefix matches
			for (map<String,String>::const_iterator it = section_descriptions.begin(); it!=section_descriptions.end(); ++it)
			{
				if (path.hasPrefix(it->first))
				{
					param_editable_->setSectionDescription(it->first, it->second);
				}
			}
			section_descriptions.clear();
		}
	
		for (Int i = 0; i < child->childCount();++i)
		{
			storeRecursive_(child->child(i),path,section_descriptions);	//whole tree recursively
		}	
	}
	  
	
	void ParamEditor::contextMenuEvent(QContextMenuEvent* event)
	{
		selected_item_ = itemAt(event->pos());
		
		QMenu menu(this);
		if(!selected_item_) // there is no item under the requested position
		{
			menu.addAction(tr("&Insert new value"), this, SLOT(insertItem()));
			menu.addAction(tr("&Insert new section"), this, SLOT(insertNode()));
		}
		else if (selected_item_->data(0,Qt::UserRole)==NODE || selected_item_->data(0,Qt::UserRole)==ITEM)
		{
			menu.addAction(tr("&Copy"), this, SLOT(copySubTree()));
			menu.addAction(tr("C&ut"), this, SLOT(cutSubTree()));
			menu.addAction(tr("&Paste"), this, SLOT(pasteSubTree()));
			menu.addSeparator();
			menu.addAction(tr("&Delete"), this, SLOT(deleteItem()));
			menu.addAction(tr("&Insert new value"), this, SLOT(insertItem()));
			menu.addAction(tr("&Insert new section"), this, SLOT(insertNode()));
		}
		if(selected_item_->data(0,Qt::UserRole)==NODE)
		{
			menu.addSeparator();
			menu.addAction(tr("&Expand"), this, SLOT(expandTree()));
			menu.addAction(tr("&Collapse"), this, SLOT(collapseTree()));
		}
		menu.exec(event->globalPos());
		selected_item_=NULL;
	}
	
	void ParamEditor::expandTree()
	{
		QTreeWidgetItem* item =currentItem();
		if(item)
		{
			expandItem(item);
		}
	}
	
	void ParamEditor::collapseTree()
	{
		QTreeWidgetItem* item =currentItem();
		if(item)
		{
			collapseItem(item);
		}
	}
	
	void ParamEditor::copySubTree()
	{
		if(!selected_item_)
		{
			selected_item_=currentItem();
			if(!selected_item_)
			{
				return;
			}
		}
		
		copied_item_=selected_item_->clone(); // the item of whom we make copies
		selected_item_=NULL; // always reset selected_item_ to NULL because we use it to check if user clicks in no item area
	}
	
	void ParamEditor::pasteSubTree()
	{
		if(!selected_item_)
		{
			selected_item_=currentItem();
			if(!selected_item_)
			{
				selected_item_=invisibleRootItem();
			}
		}
		
		if(selected_item_->data(0,Qt::UserRole)==ITEM)
		{
			selected_item_ = selected_item_->parent();
		}
		
		if(copied_item_ )
		{
			QTreeWidgetItem* new_child=copied_item_->clone();
			selected_item_->addChild(new_child);
			Int child_count=0;
			for (Int i = selected_item_->childCount()-1; i >=0;i--)
			{
				if(selected_item_->child(i)->text(0)==new_child->text(0))
				{
					if(++child_count>1)
					{
						setCurrentItem(new_child);
						editItem(new_child);
						break;
					}
				}
			}
			setModified(true);
		}
		selected_item_=NULL;
	}
	
	void ParamEditor::cutSubTree()
	{
		copySubTree();
		deleteItem();
	}
	
	void ParamEditor::setModified(bool is_modified)
	{
		if (is_modified != modified_)
		{
			modified_ = is_modified;
			emit modified(modified_);
		}
	}

	bool ParamEditor::isModified()
	{
		return modified_;
	}
		
} // namespace OpenMS
