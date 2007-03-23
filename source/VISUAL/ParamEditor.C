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

#include <QtGui/QAction>
#include <QtGui/QKeyEvent>
#include <QtGui/QMenu>


using namespace std;

namespace OpenMS
{

	ParamEditor::ParamEditor(QWidget * parent)
	  : QTreeWidget(parent),
	  	param_editable_(0),
	  	param_const_(0)
	{
		setMinimumSize(500,300);
		
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
		String path;
		if (isValid()) 
		{
			if(param_editable_!=NULL)
			{
				
				ret=true;
				QTreeWidgetItem* parent=this->invisibleRootItem();  
				param_editable_->clear();
				path+=parent->child(0)->text(0).toStdString();
				
		
				for (Int i = 0; i < parent->child(0)->childCount();++i)
				{
					storeRecursive_(parent->child(0)->child(i),path);	//whole tree recursively
				}	
			
			}
		}
		else 
		{
			std::cerr<<"\n***You must stick to the Param format!***\n";
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
		for(int i=item->childCount()-1; i>=0; i--)
		{
			deleteItemRecursive_(item->child(i));
		}
		delete item;
	}
	
	void ParamEditor::deleteAll()
	{
		QTreeWidgetItem* item=invisibleRootItem();
		deleteItemRecursive_(item->child(0));
	}
		
	void ParamEditor::insertItem()
	{
		QTreeWidgetItem* parent=currentItem();
		QTreeWidgetItem* item=NULL;
			
		if(!parent && topLevelItemCount()==0)
		{
			item=new QTreeWidgetItem(invisibleRootItem());
		}
		else if(parent && parent->text(1).size()==0 && parent->text(2)==0)
		{
			item=new QTreeWidgetItem(parent);
		}
		else
		{
			return;
		}
		item->setText(0, "");
		item->setText(1, "");
		item->setText(2, "");
		item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
		
	}
	    
	void ParamEditor::storeRecursive_(const QTreeWidgetItem* child, String path) const
	{
		path+=":"+child->text(0).toStdString();
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
		bool ret=true;
		bool ok1,ok2;
		if(parent->childCount()==0)		//   *******ITEM*******
		{
			if(parent->text(0).size()==0)
			{
				ret=false;
				
			}
			else
			{
				parent->text(0).toDouble(&ok1);
				parent->text(0).toLong(&ok2);
				if(ok1==true || ok2==true) ret=false;
				

			}
			if(parent->text(2)=="int")
			{
				parent->text(1).toLong(&ok1);
				if(!ok1) ret=false;
				if(parent->text(1).size()==0) ret=true;
				

			}
			else if(parent->text(2)=="float")
			{
				parent->text(1).toDouble(&ok1);
				if(!ok1) ret=false;
				if(parent->text(1).size()==0) ret=true;
				
			}
			else if(parent->text(2)=="string")
			{
				
				parent->text(1).toDouble(&ok1);
				parent->text(1).toLong(&ok2);
				if(ok1==true || ok2==true) ret=false;
				if(parent->text(1).size()==0) ret=true;
				
			}
			else
			{
				ret=false;
				
				
			}
		}
		else				//  *****NODE********
		{
			if(parent->text(0).size()==0 || parent->text(1).size()!=0 || parent->text(2).size()!=0)		//if name of node is empty or type and value are not empty--> wrong format
			{
				ret=false;
				
			}
			else
			{
				parent->text(0).toLong(&ok1);
				parent->text(0).toDouble(&ok2);
				//if(ok1==true || ok2==true ) ret=false;	//(if name of node can be converted to a number--> wrong format) is maybe too strict?
			}
			for (Int i = 0; i < parent->childCount();++i)
			{
				if(!(isValidRecursive_(parent->child(i)))) ret=false;	//check whole tree recursively
			}	
		}	
		return ret;
	}

	void ParamEditor::keyPressEvent(QKeyEvent* e)
	{
		if(e->key()==Qt::Key_Delete)
		{
			deleteItem();
		}
		else if(e->key()==Qt::Key_Insert)
		{
			insertItem();
		}
		else
		{
			QAbstractItemView::keyPressEvent ( e );
		}
	}

	void ParamEditor::contextMenuEvent(QContextMenuEvent* event)
	{
		QTreeWidgetItem* item = itemAt(event->pos());
		if (item  && item->text(1).size()==0 && item->text(2)==0)
		{
			QMenu menu(this);
			menu.addAction(tr("&Delete item"), this, SLOT(deleteItem()));
			menu.addAction(tr("&Insert item"), this, SLOT(insertItem()));
			menu.addAction(tr("&Expand all"), this, SLOT(expandAll()));
			menu.addAction(tr("&Collapse all"), this, SLOT(collapseAll()));			
			menu.exec(event->globalPos());
		}
		else
		{
			// there is no item under the requested position
		}
	}
	

} // namespace OpenMS
