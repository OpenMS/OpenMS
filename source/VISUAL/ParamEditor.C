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

#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/CONCEPT/Types.h>

#include <QtGui/QAction>
#include <QtGui/QKeyEvent>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>
#include <QtGui/QItemDelegate>
#include <QtCore/QModelIndex>
#include <QtGui/QComboBox>
#include <QtCore/QAbstractItemModel>
#include <QtCore/QString>
#include <QtCore/QStringList>




using namespace std;

namespace OpenMS
{
	namespace Internal
	{

		ParamEditorDelegate::ParamEditorDelegate(QObject *parent)
		     : QItemDelegate(parent)
		 {
		 }
		 
		 QWidget *ParamEditorDelegate::createEditor(QWidget *parent,
		     const QStyleOptionViewItem & option,
		     const QModelIndex & index ) const
		 {
			 QString str = index.model()->data(index, Qt::DisplayRole).toString();
			 Int id=index.model()->data(index, Qt::UserRole).toInt();
			 if(index.column()==1 && id==ParamEditor::NODE)
			 {
				 return 0;
			 }
			 else if (index.column() == 2 && id==ParamEditor::ITEM)
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
			else if(index.column()==2 && id==ParamEditor::NODE) return 0;
			 else return QItemDelegate::createEditor(parent,option,index);

		 }
		 
		 void ParamEditorDelegate::setEditorData(QWidget *editor,
						     const QModelIndex &index) const
		 {
			 if(index.column()==2)
			 {
				QString str = index.model()->data(index, Qt::DisplayRole).toString();
				 QStringList list;
				list<<"int"<<"float"<<"string";
				int pos = list.indexOf(str);
				QComboBox *combo = static_cast<QComboBox*>(editor);
				if (pos!=-1)
				{
					combo->setCurrentIndex(pos);
				}
			}
				else QItemDelegate::setEditorData(editor,index);
		 }
		 
		 void ParamEditorDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
						    const QModelIndex &index) const
		 {
			 if(index.column()==2)
			 {
				QComboBox *combo= static_cast<QComboBox*>(editor);
				QString str = combo->currentText();

				model->setData(index, str);
			 }
			 else QItemDelegate::setModelData(editor,model,index);
		 }
		 
		 void ParamEditorDelegate::updateEditorGeometry(QWidget *editor,
		     const QStyleOptionViewItem &option, const QModelIndex & index) const
		 {
			 if(index.column()==2)
			 {
				editor->setGeometry(option.rect);
			 }
			 else QItemDelegate::updateEditorGeometry(editor,option,index);
		 }
	}
 }

namespace OpenMS
{

	ParamEditor::ParamEditor(QWidget * parent)
	  : QTreeWidget(parent),
	  	param_editable_(0),
	  	param_const_(0)
	{
		setMinimumSize(500,300);
		setItemDelegate(new Internal::ParamEditorDelegate);
		setWindowTitle("ParamEditor");
		setColumnCount(3);
		QStringList list;
		list.push_back("name");
		list.push_back("value");
		list.push_back("type");
		setHeaderLabels(list);
	}
	    
	void ParamEditor::load(const Param& param)
	{
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
	//				cout << "key_wo: "<<key_without_prefix<<endl;
	//				cout << "key   : "<<key<<endl;
	//				cout << "prefix: "<<prefix<<endl;
	//				cout << "|||   : "<<key.substr(0, common)<<endl;
				//write down
			down = prefix.substr(common, prefix.size());
			if (down!="")
			{
	//				cout << "  <-  : "<<down<<endl;
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
			if (up!="")
			{
	//					cout << "  ->  : "<<up<<endl;
				while (up != "")
				{
					UInt pos = up.find(":");
					item = new QTreeWidgetItem(parent);
					item->setText(0, QString::fromStdString ( up.substr(0,pos)));
					item->setText(1, QString::fromStdString ( ""));
					item->setText(2, QString::fromStdString ( ""));
					item->setData(0,Qt::UserRole,NODE);
					item->setData(1,Qt::UserRole,NODE);
					item->setData(2,Qt::UserRole,NODE);
					if(param_editable_!=NULL)
					{
						item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
					}
					else 
					{
						item->setFlags( Qt::ItemIsEnabled );
					}
						parent=item;
						up = up.substr(pos+1,up.size());
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
	    
	bool ParamEditor::store() const
	{
		bool ret=false;
		if (isValid()) 
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
		}
		else 
		{
			QMessageBox::critical(const_cast<ParamEditor*>(this),"Error writing file!","Look at the output to stderr for details!");
		}
		return ret;
	}
	    
	bool ParamEditor::isValid() const
	{
		bool ret=true;
		QTreeWidgetItem* parent=this->invisibleRootItem();
		for (int i = 0; i < parent->childCount();++i)
			{
				if(!(isValidRecursive_(parent->child(i)))) ret=false;	//check whole tree recursively
			}	
		return ret;
	}
	    
	void ParamEditor::deleteItem()
	{
		QTreeWidgetItem* item=currentItem();
		if (!item) return;
		for(int i=item->childCount()-1; i>=0; i--)
		{
			deleteItemRecursive_(item->child(i));
		}
		delete item;
	}
		
	void ParamEditor::deleteItemRecursive_(QTreeWidgetItem* item)
	{				
		for(Int i=item->childCount()-1; i>=0; i--)
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
	}
		
	void ParamEditor::insertItem()
	{
		QTreeWidgetItem* parent=selected_item_;
		QTreeWidgetItem* item=NULL;
			
		if(!parent)
		{
			item=new QTreeWidgetItem(invisibleRootItem());
		}
		else if(parent->data(0,Qt::UserRole)==NODE)
		{
			item=new QTreeWidgetItem(parent);
		}
		else if(parent->data(0,Qt::UserRole)==ITEM)
		{
			return;
		}
		item->setText(0, "");
		item->setText(1, "");
		item->setText(2, "");
		item->setData(0,Qt::UserRole,ITEM);
		item->setData(1,Qt::UserRole,ITEM);
		item->setData(2,Qt::UserRole,ITEM);
		item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
		
	}
	    
	void ParamEditor::insertNode()
	{
		QTreeWidgetItem* parent=selected_item_;
		QTreeWidgetItem* item=NULL;
			
		if(!parent)
		{
			item=new QTreeWidgetItem(invisibleRootItem());
		}
		else if(parent->data(0,Qt::UserRole)==NODE)
		{
			item=new QTreeWidgetItem(parent);
		}
		else if(parent->data(0,Qt::UserRole)==ITEM)
		{
			return;
		}
		item->setText(0, "");
		item->setText(1, "");
		item->setText(2, "");
		item->setData(0,Qt::UserRole,NODE);
		item->setData(1,Qt::UserRole,NODE);
		item->setData(2,Qt::UserRole,NODE);
		item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
		
	}
	void ParamEditor::storeRecursive_(const QTreeWidgetItem* child, String path) const
	{
		if (path=="")
		{
			path = child->text(0).toStdString();
			
		}
		else
		{
			path = path + ":" + child->text(0).toStdString();	
		}
		
		if(child->text(2)=="float")
		{
			param_editable_->setValue(path,child->text(1).toDouble());
		}
		else if(child->text(2)=="string")
		{
			param_editable_->setValue(path, child->text(1).toStdString());
			//std::cerr<<"\n"<<path<<":  "<<child->text(1).toStdString()<<"\n";
		}
		else if(child->text(2)=="int")
		{
			param_editable_->setValue(path, child->text(1).toInt());
		}
	
		for (Int i = 0; i < child->childCount();++i)
		{
			storeRecursive_(child->child(i),path);	//whole tree recursively
		}	
	}
	   
	
	bool ParamEditor::isValidRecursive_(QTreeWidgetItem* parent) const
	{
		bool ok;
		// -- ITEM --
		if(parent->childCount()==0)		
		{
			if(parent->text(0).size()==0)
			{
				cerr << "ParamEditor: Error - Empty item name!" << endl;	
				return false;
			}
			
			String type = parent->text(2).toStdString();
			if(type=="int")
			{
				parent->text(1).toLong(&ok);
				if(!ok)
				{
					cerr << "ParamEditor: Error - Invalid 'int' value '" << parent->text(1).toStdString() << "'!" << endl;	
					return false;
				}
			}
			else if(type=="float")
			{
				parent->text(1).toDouble(&ok);
				if(!ok)
				{
					cerr << "ParamEditor: Error - Invalid 'float' value '" << parent->text(1).toStdString() << "'!" << endl;	
					return false;
				}
			}
			else if(type=="string")
			{
				//nothing to check here, but not removable because of else case
			}
			else
			{
				cerr << "ParamEditor: Error - Unknown type '" << type << "'!" << endl;	
				return false;
			}
		}
		// -- NODE --
		else 
		{
			//name of node is empty
			if(parent->text(0).size()==0)		
			{
				cerr << "ParamEditor: Error - Empty subsection!" << endl;	
				return false;
			}
			//check whole tree recursively
			ok = true;
			for (Int i = 0; i < parent->childCount();++i)
			{
				if(isValidRecursive_(parent->child(i))==false)
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

	void ParamEditor::keyPressEvent(QKeyEvent* e)
	{
		if(e->key()==Qt::Key_Delete)
		{
			selected_item_=currentItem();
			deleteItem();
		}
		else if(e->key()==Qt::Key_Insert)
		{
			selected_item_=currentItem();
			insertItem();
		}
		else
		{
			QAbstractItemView::keyPressEvent ( e );
		}
	}

	void ParamEditor::contextMenuEvent(QContextMenuEvent* event)
	{
		selected_item_ = itemAt(event->pos());
		if(!selected_item_)
		{
			// there is no item under the requested position
			QMenu menu(this);
			menu.addAction(tr("&Insert item"), this, SLOT(insertItem()));
			menu.addAction(tr("&Insert node"), this, SLOT(insertNode()));
			menu.exec(event->globalPos());
		}
		else if (selected_item_->data(0,Qt::UserRole)==NODE)
		{
			QMenu menu(this);
			menu.addAction(tr("&Delete item"), this, SLOT(deleteItem()));
			menu.addAction(tr("&Insert item"), this, SLOT(insertItem()));
			menu.addAction(tr("&Insert node"), this, SLOT(insertNode()));
			menu.addAction(tr("&Expand subtree"), this, SLOT(expandTree()));
			menu.addAction(tr("&Collapse subtree"), this, SLOT(collapseTree()));			
			menu.exec(event->globalPos());
		}
		else if(selected_item_->data(0,Qt::UserRole)==ITEM)
		{
			QMenu menu(this);
			menu.addAction(tr("&Delete item"), this, SLOT(deleteItem()));
			menu.exec(event->globalPos());
		}
		
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
	
	

} // namespace OpenMS
