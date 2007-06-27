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
#include <OpenMS/FORMAT/Param.h>

#include <QtGui/QAction>
#include <QtGui/QKeyEvent>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>
#include <QtGui/QItemDelegate>
#include <QtCore/QModelIndex>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>
#include <QtCore/QAbstractItemModel>
#include <QtCore/QString>
#include <QtCore/QStringList>
#include <QtGui/QRegExpValidator>
#include <QtCore/QRegExp>
#include <QtGui/QValidator>
#include <QtGui/QColor>
#include <QtGui/QPainter>
#include <QtGui/QBrush>
#include <QtGui/QShortcut>


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
			QString str = index.model()->data(index, Qt::DisplayRole).toString();
			Int id=index.model()->data(index, Qt::UserRole).toInt();
			if(index.column()==0)	// if it is the first column give me a QLineEdit with a regex validator that allows only words
			{
				QRegExp rx("\\S+");
				QValidator *validator = new QRegExpValidator(rx, parent);
				QLineEdit *editor = new QLineEdit(parent);
				editor->setValidator(validator);
				return editor;
			}
			else if (index.column() == 2 && id==ParamEditor::ITEM) // if we have an item and we are in the third column give me a combobox to choose the type 
			{
				QComboBox *editor = new QComboBox(parent);
				QStringList list;
				list<<"int"<<"float"<<"string";
				editor->addItems(list);
				int pos =list.indexOf(str);
				if (pos!=-1)
				{
					editor->setCurrentIndex(pos);
				}
				return editor;
			}
			else if(index.column()==1 && id==ParamEditor::ITEM)	// if we are in the second column and it's an item give me a QLineEdit
			{
				QLineEdit *editor = new QLineEdit(parent);
				return editor;
			}

			return 0;
		}
		 
		void ParamEditorDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
		{
			if(index.column()==2)
			{
				QString str = index.model()->data(index, Qt::DisplayRole).toString();
				QStringList list;
				list<<"int"<<"float"<<"string";
				int pos = list.indexOf(str);	// look for the string that is current displayed type
				QComboBox *combo = static_cast<QComboBox*>(editor);
				if (pos!=-1)
				{
					combo->setCurrentIndex(pos);	//set the current type in the combobox
				}
			}
			else if(index.column()==1)
			{
				QString str = index.model()->data(index, Qt::DisplayRole).toString();	// get the value from the model
				QLineEdit *edit = static_cast<QLineEdit*>(editor);	// cast the editor to your actual widget
				edit->setText(str);	// set it to the current value from the model
			}
			else if(index.column()==0)
			{
				QString str = index.model()->data(index, Qt::DisplayRole).toString();
				QLineEdit *edit = static_cast<QLineEdit*>(editor);
				edit->setText(str);
		 	}
		}
		 
		void ParamEditorDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
		{
			QVariant present_value = index.model()->data(index, Qt::DisplayRole);
			
			if(index.column()==2)
			{
				
				QComboBox *combo= static_cast<QComboBox*>(editor);
				QVariant new_value(combo->currentText());
				
				if(new_value!=present_value)	// if data of a cell gets a new value from the editor widget change color to yellow and emit the modified signal
				{
					model->setData(index,QBrush(Qt::yellow),Qt::BackgroundRole);
					emit modified(true);
				}
				model->setData(index, new_value);
			}
			else if(index.column()==1)
			{
				QLineEdit *edit= static_cast<QLineEdit*>(editor);
				QVariant new_value(edit->text());
				
				if(new_value!=present_value)
				{
					model->setData(index,QBrush(Qt::yellow),Qt::BackgroundRole);
					emit modified(true);
				}
				model->setData(index, new_value);
			}
			else if(index.column()==0)
			{
				QLineEdit *edit= static_cast<QLineEdit*>(editor);
				QVariant new_value(edit->text());
			
				if(new_value!=present_value)
				{
					model->setData(index,QBrush(Qt::yellow),Qt::BackgroundRole);
					emit modified(true);
				}
				model->setData(index, new_value);
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
			modified_(false),
			modificationsCount_(0),
			is_name_empty_(false)
	{
		setMinimumSize(500,300);
		setItemDelegate(new Internal::ParamEditorDelegate);	// the delegate from above is set
		setWindowTitle("ParamEditor");
		setColumnCount(3);
		connect(itemDelegate(),SIGNAL(modified(bool)),this,SLOT(setModified(bool)));	// the modified signal from the delegate is connected with setModified that operates a counter variable to watch over changes
		connect(this, SIGNAL(currentItemChanged ( QTreeWidgetItem *, QTreeWidgetItem*)),this,SLOT(editChanged(QTreeWidgetItem*, QTreeWidgetItem*))); 
		QStringList list;
		list.push_back("name");
		list.push_back("value");
		list.push_back("type");
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
	
	void ParamEditor::editChanged(QTreeWidgetItem* current, QTreeWidgetItem* previous)
	{
		if(previous)	// if we have a previous item, if not previous is NULL and we have a crash
		{
			QTreeWidgetItem* parent=previous->parent();	

			if(!parent)
			{
				parent=invisibleRootItem();
			}
			Int child_count=0;
			for (Int i = parent->childCount()-1; i >=0;i--)	// we make sure that all childs of a section have a unique name if not editItem() when left for a new item
			{
				if(parent->child(i)->text(0)==previous->text(0))
				{
					++child_count;
					if(child_count>1 && current->text(0)!=previous->text(0))
					{
						setCurrentItem(previous);
						editItem(previous);
						break;
					}
					else if(child_count>1 && current->text(0)==previous->text(0))
					{
						editItem(current);
						break;
					}
				}
			}
			if(previous->text(0).isEmpty()) // we make sure that previous item has a name if not editItem()
			{	
				editItem(previous);
				setCurrentItem(previous);
				is_name_empty_=true;
			}
			else
			{
				is_name_empty_=false;
			}
		}
	}

	void ParamEditor::load(const Param& param)
	{
		clear();
		
		string up, down ,key, key_without_prefix, new_prefix ,type, prefix = "";
		UInt common;//, level=1;
		QTreeWidgetItem* parent=this->invisibleRootItem();
		param_const_=&param;
		QTreeWidgetItem* item = NULL;	
		for(Param::ConstIterator it=param.begin();it!=param.end();++it)
		{
			//init variables
			key = it->first;
			if (key.find(":")==string::npos)
			{
				key_without_prefix = key;
				new_prefix = "";
			}
			else
			{
				key_without_prefix = key.substr(key.rfind(":")+1,key.size());
				new_prefix = key.substr(0,key.rfind(":")+1);
			} 
			
			//common prefix
			common=0;
			for (UInt i=0;i<min(key.size(),prefix.size());++i)	
			{											
				if (prefix[i]!=key[i])
				{
					break;
				}
				if (prefix[i]==':')							
				{
					common = i+1;
				}
			}
			//cout << "key_wo: "<<key_without_prefix<<endl;
			//cout << "key   : "<<key<<endl;
			//cout << "prefix: "<<prefix<<endl;
			//cout << "|||   : "<<key.substr(0, common)<<endl;
			
			//write down
			down = prefix.substr(common, prefix.size());
			if (down!="")
			{
				//cout << "  <-  : "<<down<<endl;
				for (UInt i = 0; i < down.size();++i)	
				{
					if (down[i]==':')
					{
						parent=parent->parent();
						if(parent==NULL) parent=invisibleRootItem();
					}
				}	
			}
				
			//write up
			up = key.substr(common, key.size()-common-key_without_prefix.size());
			String nodepath = key.substr(0,common);
			if (!nodepath.empty())
			{
				nodepath = nodepath.substr(0,-1);
			}
			if (up!="")
			{
				//cout << "  ->  : "<<up<<endl;
				while (up != "")
				{
					String nodename = up.substr(0,up.find(":"));
					item = new QTreeWidgetItem(parent);
					item->setText(0, QString::fromStdString ( nodename));
					item->setText(1, QString::fromStdString ( ""));
					item->setText(2, QString::fromStdString ( ""));
					item->setData(0,Qt::UserRole,NODE);
					item->setData(1,Qt::UserRole,NODE);
					item->setData(2,Qt::UserRole,NODE);
					
					//description
					if (nodepath.empty())
					{
						nodepath = nodename;
					}
					else
					{
						nodepath = nodepath + ":" + nodename;
					}
					//cout << "NODE: '" << nodepath << "': " << param.getDescription(nodepath) << endl;
					String description = param.getDescription(nodepath);
					description.substitute("\n","<BR>");
					//cout << description << endl << endl;
					item->setToolTip(0,description.toQString());
					
					//flags
					if(param_editable_!=NULL)
					{
						item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
					}
					else 
					{
						item->setFlags( Qt::ItemIsEnabled );
					}
					
					//update loop variables
					parent=item;
					up = up.substr(nodename.size()+1,up.size());
				}				
			}
				
			if (it->second.valueType()==DataValue::INTVALUE || it->second.valueType()==DataValue::LONVALUE || it->second.valueType()==DataValue::SHOVALUE  )
			{
				type = "int";
			}
			
			if (it->second.valueType()==DataValue::FLOVALUE || it->second.valueType()==DataValue::DOUVALUE )
			{
				type = "float";
			}
			
			if (it->second.valueType()==DataValue::STRVALUE )
			{
				type = "string";
			}
				
			if(it->second.valueType()!=DataValue::EMPTYVALUE)
			{
				item = new QTreeWidgetItem(parent);
				item->setText(0, QString::fromStdString ( key_without_prefix));
				item->setText(1, QString::fromStdString ( it->second.toString()));
				item->setText(2, QString::fromStdString ( type));
				String description = param.getDescription(it->first);
				description.substitute("\n","<BR>");
				item->setToolTip(0,description.toQString());
				//cout << description << endl << endl;
				//cout << "ITEM: '" << key_without_prefix << "': " << param.getDescription(it->first) << endl;
				item->setData(0,Qt::UserRole,ITEM);
				item->setData(1,Qt::UserRole,ITEM);
				item->setData(2,Qt::UserRole,ITEM);
				if(param_editable_!=NULL)
				{
					item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
				}
				else 
				{
					item->setFlags( Qt::ItemIsEnabled );
				}
			}
				
				//set new prefix
			prefix = new_prefix;
		}
			
			//close remaining prefix tags
			down = prefix;
		if (down!="")
		{
	//			cout << "  <-  : "<<down<<endl;
			for (UInt i = 0; i < down.size();++i)
			{
				if (down[i]==':')
				{
					parent=parent->parent();
					if(parent==NULL) parent=invisibleRootItem();
				}
			}	
		}
		expandAll();
		
		resizeColumnToContents(0);
		resizeColumnToContents(1);
		resizeColumnToContents(2);
		
	}

	void ParamEditor::loadEditable(Param& param)
	{
		param_editable_=&param;
		load(param);
	}
	    
	bool ParamEditor::store()
	{
		bool ret=false;
		QStringList list;
		if (isValid(list)) 
		{
			if(param_editable_!=NULL)
			{
				ret=true;
				QTreeWidgetItem* parent=this->invisibleRootItem();
				param_editable_->clear();
			
				for (Int i = 0; i < parent->childCount();++i)
				{
					storeRecursive_(parent->child(i),"");	//whole tree recursively
				}	
			}
			modificationsCount_=1;
			setModified(false);
		}
		else 
		{
			QMessageBox::critical(this,"ParamEditor: Error writing file!",list.join("\n"));
		}
		return ret;
	}
	    
	bool ParamEditor::isValid(QStringList& list) const
	{
		bool ret=true;
		QTreeWidgetItem* parent=this->invisibleRootItem();
		for (int i = 0; i < parent->childCount();++i)
			{
				if(!(isValidRecursive_(parent->child(i),list))) ret=false;	//check whole tree recursively
			}	
			
		return ret;
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
		item->setText(0, "");
		item->setText(1, "");
		item->setText(2, "string");
		item->setData(0,Qt::UserRole,ITEM);
		item->setData(1,Qt::UserRole,ITEM);
		item->setData(2,Qt::UserRole,ITEM);
		item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
		item->setBackground (0,Qt::yellow);
		item->setBackground (1,Qt::yellow);
		item->setBackground (2,Qt::yellow);
		setCurrentItem(item);	// we set the new item to be current to prevent an empty name by editChanged
		editItem(item);			// we edit the item because it is empty
		setModified(true);		// and we tell everyone it's been modified
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
		item->setText(0, "");
		item->setText(1, "");
		item->setText(2, "");
		item->setData(0,Qt::UserRole,NODE);
		item->setData(1,Qt::UserRole,NODE);
		item->setData(2,Qt::UserRole,NODE);
		item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
		item->setBackground (0,Qt::yellow);
		item->setBackground (1,Qt::yellow);
		item->setBackground (2,Qt::yellow);
		setCurrentItem(item);
		editItem(item);
		setModified(true);
	}
	
	void ParamEditor::storeRecursive_(QTreeWidgetItem* child, String path)
	{
		child->setData ( 0, Qt::BackgroundRole, QBrush(Qt::white));
		child->setData ( 1, Qt::BackgroundRole, QBrush(Qt::white));
		child->setData ( 2, Qt::BackgroundRole, QBrush(Qt::white));
		
		if (path=="")
		{
			path = child->text(0).toStdString();
			
		}
		else
		{
			path = path + ":" + child->text(0).toStdString();	
		}
		
		//write back the description
		String description = child->toolTip(0);
		description.substitute("<BR>","\n");
		
		if(child->text(2)=="float")
		{
			param_editable_->setValue(path,child->text(1).toDouble(),description);
		}
		else if(child->text(2)=="string")
		{
			param_editable_->setValue(path, child->text(1).toStdString(),description);
			//std::cerr<<"\n"<<path<<":  "<<child->text(1).toStdString()<<"\n";
		}
		else if(child->text(2)=="int")
		{
			param_editable_->setValue(path, child->text(1).toInt(),description);
		}
	
		for (Int i = 0; i < child->childCount();++i)
		{
			storeRecursive_(child->child(i),path);	//whole tree recursively
		}	
	}
	   
	
	bool ParamEditor::isValidRecursive_(QTreeWidgetItem* parent,QStringList& list) const
	{
		bool ok;
		// -- ITEM --
		if(parent->data(0,Qt::UserRole)==ITEM)		
		{
			if(parent->text(0).size()==0)
			{
				list<< "Error - Empty item name!";	
				return false;
			}
			
			String type = parent->text(2).toStdString();
			if(type=="int")	// we check if the type of an item is okay with the actual value, we do this by QString's conversion methods
			{
				parent->text(1).toLong(&ok);
				if(!ok)
				{
					list<<QString("Error - Invalid 'int' value '%1'!").arg(parent->text(1));	
					return false;
				}
			}
			else if(type=="float")
			{
				parent->text(1).toDouble(&ok);
				if(!ok)
				{
					list<< QString("Error - Invalid 'float' value '%1'!").arg(parent->text(1));	
					return false;
				}
			}
			else if(type=="string")
			{
				//nothing to check here, but not removable because of else case
			}
			else
			{
				list<< QString("Error - Unknown type '%1'!").arg(type.c_str());	
				return false;
			}
		}
		// -- NODE --
		else 
		{
			//name of node is empty
			if(parent->text(0).size()==0)		
			{
				list<< "Error - Empty subsection!";	
				return false;
			}
			//check whole tree recursively
			ok = true;
			for (Int i = 0; i < parent->childCount();++i)
			{
				if(isValidRecursive_(parent->child(i), list)==false)
				{
					ok = false;
				}
			}
			if (!ok)
			{
				return false;
			}
		}	
		return true;
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
			if(is_modified)	// here we increment or decement changes in the model that are signaled by the signal modified to see when the data is back in original state
			{
				++modificationsCount_;
				//cerr<<"\n"<<modificationsCount_<<"\n";
			}
			else 
			{
				--modificationsCount_;
				//cerr<<"\n"<<modificationsCount_<<"\n";

			}
	
		if(modificationsCount_==0)
		{
			modified_=false;	// == "original state"
		}
		else
		{
			modified_=true;
		}
		emit modified(modified_); // ParamEditor itself emits a modified signal for other widgets to check on the changes of ParamEditor
	}

	bool ParamEditor::isModified()
	{
		return modified_;
	}
	
	bool ParamEditor::isNameEmpty()
	{
		return is_name_empty_;
	}
	
} // namespace OpenMS
