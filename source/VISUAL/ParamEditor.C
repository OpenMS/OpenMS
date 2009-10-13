// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/VISUAL/ListEditor.h>

#include <QtGui/QMessageBox>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>
#include <QtGui/QShortcut>
#include <QtGui/QMenu>
#include <QtGui/QItemSelection>
#include <QtCore/QStringList>
#include <QtGui/QLabel>
#include <QtGui/QFileDialog>
 
#include <stack>
#include <limits>
#include <sstream>

using namespace std;

/*

Description of the data stored in the itmes:

            | Column 0  | Column 1    | Column 2 | Column 3         |
---------------------------------------------------------------------  
DisplayRole | name      | value       | type     | restr. (display) |
UserRole    | NODE/ITEM | description | restr.   |                  |


*/

namespace OpenMS
{
	namespace Internal
	{
		ParamEditorDelegate::ParamEditorDelegate(QObject* parent)
			: QItemDelegate(parent)
		{
		}
		
		QWidget *ParamEditorDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem& , const QModelIndex& index) const
		{
			Int type = index.sibling(index.row(),0).data(Qt::UserRole).toInt();
			if(index.column()==1 && type!=ParamEditor::NODE)
			{
				QString dtype = index.sibling(index.row(),2).data(Qt::DisplayRole).toString();
				QString restrictions = index.sibling(index.row(),2).data(Qt::UserRole).toString();
				if (dtype=="string" && restrictions!="") //Drop-down list for enums
				{
					QComboBox* editor = new QComboBox(parent);
					QStringList list;
					list.append("");
					list += restrictions.split(",");
					editor->addItems(list);
					connect(editor,SIGNAL(activated(int)),this, SLOT(commitAndCloseComboBox_()));
					return editor;
				}
				else if(dtype == "output file")
				{
					QLineEdit* editor = new QLineEdit(parent);
					//editor->setReadOnly(true);
					QString str = index.sibling(index.row(),0).data(Qt::DisplayRole).toString();
					fileName_ = QFileDialog::getSaveFileName(editor,tr("Output File"),str);
					return editor;
				}
				else if(dtype == "input file")
				{
					QLineEdit* editor = new QLineEdit(parent);
					//editor->setReadOnly(true);
					QString str = index.sibling(index.row(),0).data(Qt::DisplayRole).toString();
					fileName_ = QFileDialog::getOpenFileName(editor,tr("Input File"),str);
					return editor;
				}
				else if(dtype =="string list" || dtype =="int list" || dtype=="double list" || dtype =="input file list" || dtype =="output file list" ) // for lists
				{
					QString str = "<"+ index.sibling(index.row(),0).data(Qt::DisplayRole).toString() + "> " +"(<" + dtype + ">)";
					ListEditor *editor = new ListEditor(0,str);
					editor->setTypeName(index.sibling(index.row(),0).data(Qt::DisplayRole).toString());
					editor->setModal(true);
					connect(editor,SIGNAL(accepted()),this, SLOT(commitAndCloseListEditor_()));
					connect(editor,SIGNAL(rejected()),this, SLOT(closeListEditor_()));
					return editor;
				}						
				else //LineEditor for rest
				{
					QLineEdit *editor = new QLineEdit(parent);
					editor->setFocusPolicy(Qt::StrongFocus);
					return editor;
				}
			}
			return 0;
		}
		
		void ParamEditorDelegate::setEditorData(QWidget* editor, const QModelIndex& index) const
		{
			QString str = index.data(Qt::DisplayRole).toString();

			if(index.column()==1)
			{
				if(qobject_cast<QComboBox*>(editor)) //Drop-down list for enums
				{
					int index = static_cast<QComboBox*>(editor)->findText(str);
					if (index==-1)
					{
						index = 0;
					}
					static_cast<QComboBox*>(editor)->setCurrentIndex(index);
				}
				else if(qobject_cast<QLineEdit*>(editor))// LineEdit for other values
				{
					QString dtype = index.sibling(index.row(),2).data(Qt::DisplayRole).toString();
					if(dtype == "output file" || dtype == "input file") /// for output/input file
					{
						if(!fileName_.isNull())
						{
							static_cast<QLineEdit*>(editor)->setText(fileName_);
						}
					}
					else
					{
						if(str == "" && (dtype == "int" || dtype == "float"))
						{
							static_cast<QLineEdit*>(editor)->setText("0");
						}
						else
						{
							static_cast<QLineEdit*>(editor)->setText(str);
						}
					}
				}
				else		// ListEditor for lists 
				{
					String list;
					list = str.mid(1,str.length()-2);
					QString type = index.sibling(index.row(),2).data(Qt::DisplayRole).toString();
					StringList rlist = StringList::create(list);
					for(UInt i = 0; i < rlist.size(); ++i)
					{
						rlist[i]  = rlist[i].trim();
					}
					if(type == "int list")
					{
						static_cast<ListEditor*>(editor)->setList(rlist,ListEditor::INT);
					}
					else if(type == "double list")
					{
						static_cast<ListEditor*>(editor)->setList(rlist,ListEditor::FLOAT);
					}
					else if(type == "string list")
					{
						static_cast<ListEditor*>(editor)->setList(rlist,ListEditor::STRING);
					}
					else if(type == "input file list")
					{
						static_cast<ListEditor*>(editor)->setList(rlist,ListEditor::INPUT_FILE);
					}
					else if(type == "output file list")
					{
						static_cast<ListEditor*>(editor)->setList(rlist,ListEditor::OUTPUT_FILE);
					}
					static_cast<ListEditor*>(editor)->setListRestrictions(index.sibling(index.row(),2).data(Qt::UserRole).toString());
				}
			}
		}
		
		void ParamEditorDelegate::setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const
		{
			QVariant present_value = index.data(Qt::DisplayRole);
			QVariant new_value;
			StringList list;
			bool new_list = false;
			if(index.column()==1)
			{
				//extract new value
				if(qobject_cast<QComboBox*>(editor)) //Drop-down list for enums
				{
					new_value = QVariant(static_cast<QComboBox*>(editor)->currentText());
				}
				else if(qobject_cast<QLineEdit*>(editor))
				{
					QString dtype = index.sibling(index.row(),2).data(Qt::DisplayRole).toString();
					if(dtype =="output file" || dtype == "input file")// input/outut file 
					{

						new_value = QVariant(static_cast<QLineEdit*>(editor)->text());
						fileName_ = "\0";
					}
					else if(static_cast<QLineEdit*>(editor)->text() == "" &&(dtype == "int" || dtype == "float"))//numeric
					{
						new_value = QVariant("0");
					}
					else
					{
						new_value = QVariant(static_cast<QLineEdit*>(editor)->text());
					}
				}
				else
				{
					list = static_cast<ListEditor*>(editor)->getList();
					for(UInt i = 1; i < list.size(); ++i)
					{
						list[i] = "\n" + list[i] ;
					}
					new_list = true;
				}
				//check if it matches the restrictions or is empty
				if (new_value.toString()!="")
				{
					QString type = index.sibling(index.row(),2).data(Qt::DisplayRole).toString();
					bool restrictions_met = true;
					String restrictions = index.sibling(index.row(),2).data(Qt::UserRole).toString();
					if (type=="int") //check if valid integer
					{
						bool ok(true);
						new_value.toString().toLong(&ok);
						if (!ok)
						{
							QMessageBox::warning(0,"Invalid value",QString("Cannot convert '%1' to integer number!").arg(new_value.toString()) );
							new_value = present_value;
						}
						//restrictions
						vector<String> parts;
						if (restrictions.split(' ',parts))
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
					else if (type=="float") //check if valid float
					{
						bool ok(true);
						new_value.toString().toDouble(&ok);
						if (!ok)
						{
							QMessageBox::warning(0,"Invalid value",QString("Cannot convert '%1' to floating point number!").arg(new_value.toString()) );
							new_value = present_value;
						}
						//restrictions
						vector<String> parts;
						if (restrictions.split(' ',parts))
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
					if(!restrictions_met)
					{
						QMessageBox::warning(0,"Invalid value",QString("Value restrictions not met: %1").arg(index.sibling(index.row(),3).data(Qt::DisplayRole).toString()) );
						new_value = present_value;
					}
				}
			}
			if(new_list)
			{
				stringstream ss;
				ss << list;
				QVariant new_value;
				new_value = QVariant(QString::fromStdString(ss.str()));
				model->setData(index, new_value);
				model->setData(index,QBrush(Qt::yellow),Qt::BackgroundRole);
				emit modified(true);
			}				
			else
			{
				//check if modified
				if(new_value!=present_value)
				{
					model->setData(index, new_value);
					model->setData(index,QBrush(Qt::yellow),Qt::BackgroundRole);
					emit modified(true);
				}
			}
		}
		 
		void ParamEditorDelegate::updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex& ) const
		{
				editor->setGeometry(option.rect);
		}

		bool ParamEditorDelegate::exists_(QString name, QModelIndex index) const
		{
			UInt current_index = 0;
			while(index.parent().child(current_index,0).isValid())
			{
				if (
						current_index != (UInt)(index.row())
						&& 
						index.parent().child(current_index,0).data(Qt::DisplayRole).toString()==name
						&&
						(
						(index.data(Qt::UserRole).toInt()==0 && index.parent().child(current_index,0).data(Qt::UserRole).toInt()==0)
						||
						(index.data(Qt::UserRole).toInt()!=0 && index.parent().child(current_index,0).data(Qt::UserRole).toInt()!=0)
						)
					 )
				{
					return true;
				}
				++current_index;
			}
			return false;
		}
		
		void ParamEditorDelegate::commitAndCloseListEditor_()
		{
			ListEditor* editor = qobject_cast<ListEditor*>(sender());
			emit commitData(editor);
			emit closeEditor(editor);
		}

		void ParamEditorDelegate::commitAndCloseComboBox_()
		{
			QComboBox* editor = qobject_cast<QComboBox*>(sender());
			emit commitData(editor);
			emit closeEditor(editor);
		}
		
		void ParamEditorDelegate::closeListEditor_()
		{
			ListEditor* editor = qobject_cast<ListEditor*>(sender());
			emit closeEditor(editor);
		}

		///////////////////ParamTree/////////////////////////////////
	
		ParamTree::ParamTree(QWidget* parent)
			: QTreeWidget(parent)
		{
		}
	
		void ParamTree::selectionChanged(const QItemSelection& s, const QItemSelection&)
		{
			if (!s.empty())
			{
				emit selected(s.indexes().first());
			}
		}

		bool ParamTree::edit(const QModelIndex& index, EditTrigger trigger, QEvent* event)
	 	{
			if (trigger == QAbstractItemView::EditKeyPressed)
			{
				return QAbstractItemView::edit(index.sibling(index.row(),1), trigger, event);
			}
			return QAbstractItemView::edit(index, trigger, event);
		}

	}

	///////////////////ParamEditor/////////////////////////////////

	ParamEditor::ParamEditor(QWidget* parent)
	  : QWidget(parent),
	  	param_(0),
			modified_(false),
			advanced_mode_(false)
	{
		setupUi(this);
		tree_ = new Internal::ParamTree(this);
		tree_->setMinimumSize(450,200);
		tree_->setAllColumnsShowFocus(true);
		tree_->setColumnCount(4);
		QStringList list;
		list << "name" << "value" << "type" << "restrictions";
		tree_->setHeaderLabels(list);
		dynamic_cast<QVBoxLayout*>(layout())->insertWidget(0,tree_,1);
		tree_->setItemDelegate(new Internal::ParamEditorDelegate(tree_));	// the delegate from above is set
		connect(tree_->itemDelegate(),SIGNAL(modified(bool)),this,SLOT(setModified(bool)));
		connect(advanced_,SIGNAL(toggled(bool)),this,SLOT(toggleAdvancedMode(bool)));
		connect(tree_,SIGNAL(selected(const QModelIndex&)),this,SLOT(showDocumentation(const QModelIndex&)));
	}
	
	void ParamEditor::showDocumentation(const QModelIndex& index)
	{
		doc_->setText(index.sibling(index.row(),1).data(Qt::UserRole).toString());
	}

	void ParamEditor::load(Param& param)
	{
		param_= &param;
		
		tree_->clear();
		
		QTreeWidgetItem* parent=tree_->invisibleRootItem();
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
					//item->setTextAlignment(0,Qt::AlignTop);
					//item->setTextAlignment(1,Qt::AlignTop);
					//item->setTextAlignment(2,Qt::AlignTop);
					//item->setTextAlignment(3,Qt::AlignTop);
					//name
					item->setText(0, it2->name.toQString());
					//description
					item->setData(1,Qt::UserRole,it2->description.toQString());
					//role
					item->setData(0,Qt::UserRole,NODE);
					//flags
					if(param_!=NULL)
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
					if(parent==NULL) parent=tree_->invisibleRootItem();
				}
			}
			
			//********handle item********
			item = new QTreeWidgetItem(parent);
			//item->setTextAlignment(0,Qt::AlignTop);
			//item->setTextAlignment(1,Qt::AlignTop);
			//item->setTextAlignment(2,Qt::AlignTop);
			//item->setTextAlignment(3,Qt::AlignTop);
			if (it->tags.count("advanced"))
			{
				item->setData(0,Qt::UserRole,ADVANCED_ITEM);
			}
			else //advanced parameter
			{
				item->setData(0,Qt::UserRole,NORMAL_ITEM);					
			}
			//name
			item->setText(0, it->name.toQString());
			//value
			if(it->value.valueType() == DataValue::STRING_LIST)
			{
				StringList string_list = it->value;
				String list_string = String("[")+string_list.concatenate(",\n")+"]";
				item->setText(1, list_string.toQString());
			}
			else if(it->value.valueType() == DataValue::INT_LIST)
			{
				IntList list = it->value;
				StringList string_list;
				for (Size i=0; i<list.size();++i)
				{
					string_list.push_back(list[i]);
				}
				String list_string = String("[")+string_list.concatenate(",\n")+"]";
				item->setText(1, list_string.toQString());
			}
			else if(it->value.valueType() == DataValue::DOUBLE_LIST)
			{
				DoubleList list = it->value;
				StringList string_list;
				for (Size i=0; i<list.size();++i)
				{
					string_list.push_back(list[i]);
				}
				String list_string = String("[")+string_list.concatenate(",\n")+"]";
				item->setText(1, list_string.toQString());
			}
			else
			{
				item->setText(1, String(it->value).toQString());
			}
			//type
			switch(it->value.valueType())
			{
				case DataValue::INT_VALUE:
					item->setText(2, "int");
					break;
				case DataValue::DOUBLE_VALUE:
					item->setText(2, "float");
					break;
				case DataValue::STRING_VALUE:
					if(it->tags.count("input file"))
					{
						item->setText(2,"input file");
					}
					else if(it->tags.count("output file"))
					{
						item->setText(2,"output file");
					}
					else
					{
						item->setText(2,"string");
					}
					break;
				case DataValue::STRING_LIST:
					if(it->tags.count("input file"))
					{
						item->setText(2,"input file list");
					}
					else if(it->tags.count("output file"))
					{
						item->setText(2,"output file list");
					}
					else
					{
						item->setText(2,"string list");
					}
					break;
				case DataValue::INT_LIST:
					item->setText(2,"int list");
					break;
				case DataValue::DOUBLE_LIST:
					item->setText(2,"double list");
					break;
				default:
					break;
			};
			//restrictions (displayed and internal for easier parsing)
			switch(it->value.valueType())
			{
				case DataValue::INT_VALUE:
				case DataValue::INT_LIST:
					{
						String drest="", irest="";
						bool min_set = (it->min_int!=-numeric_limits<Int>::max());
						bool max_set = (it->max_int!=numeric_limits<Int>::max());
						if (max_set || min_set)
						{
							if (min_set)
							{
								drest += String("min: ") + it->min_int;
								irest += it->min_int;
							}
							irest += " ";
							if (max_set)
							{
								if (min_set && max_set) drest += " ";
								drest += String("max: ") + it->max_int;
								irest += it->max_int;
							}
							item->setText(3, drest.toQString());
						}
						item->setData(2,Qt::UserRole,irest.toQString());
					}
					break;
				case DataValue::DOUBLE_VALUE:
				case DataValue::DOUBLE_LIST:
					{
						String drest="", irest="";
						bool min_set = (it->min_float!=-numeric_limits<DoubleReal>::max());
						bool max_set = (it->max_float!=numeric_limits<DoubleReal>::max());
						if (max_set || min_set)
						{
							if (min_set)
							{
								drest += String("min: ") + it->min_float;
								irest += it->min_float;
							}
							irest += " ";
							if (max_set)
							{
								if (min_set && max_set) drest += " ";
								drest += String("max: ") + it->max_float;
								irest += it->max_float;
							}
							item->setText(3, drest.toQString());
						}
						item->setData(2,Qt::UserRole,irest.toQString());
					}
					break;
				case DataValue::STRING_VALUE:
				case DataValue::STRING_LIST:	
					{
						String irest="";
						if (it->valid_strings.size()!=0)
						{
							irest.concatenate(it->valid_strings.begin(),it->valid_strings.end(),",");
						}
						if (irest!="")
						{
							item->setText(3, irest.toQString());
						}
						item->setData(2,Qt::UserRole,irest.toQString());
					}
					break;
				default:
					break;
			};

			//description
			item->setData(1,Qt::UserRole,it->description.toQString());
			//flags
			if(param_!=NULL)
			{
				item->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
			}
			else 
			{
				item->setFlags( Qt::ItemIsEnabled );
			}
		}
		
		tree_->expandAll();
		toggleAdvancedMode(advanced_mode_);
		
		tree_->resizeColumnToContents(0);
		tree_->resizeColumnToContents(1);
		tree_->resizeColumnToContents(2);
		tree_->resizeColumnToContents(3);
	}
	    
	void ParamEditor::store()
	{
		if(param_!=NULL)
		{
			QTreeWidgetItem* parent=tree_->invisibleRootItem();
			param_->clear();
		
			for (Int i = 0; i < parent->childCount();++i)
			{
				map<String,String> section_descriptions;
				storeRecursive_(parent->child(i),"", section_descriptions);	//whole tree recursively
			}	
		}
			
		setModified(false);
	}
	
	void ParamEditor::clear()
	{
		tree_->clear();
	}
	
	void ParamEditor::storeRecursive_(QTreeWidgetItem* child, String path,map<String,String>& section_descriptions)
	{
		child->setData ( 1, Qt::BackgroundRole, QBrush(Qt::white));
		
		if (path=="")
		{
			path = child->text(0).toStdString();
		}
		else
		{
			path += String(":") + String(child->text(0).toStdString());	
		}
		
		String description = child->data(1, Qt::UserRole).toString();
		
		if(child->text(2)=="") // node
		{
			if (description != "")
			{
				section_descriptions.insert(make_pair(path,description));
			}
		}
		else //item + section descriptions
		{
			StringList tags;
			if (child->data(0, Qt::UserRole)==ADVANCED_ITEM)
			{
				tags.push_back("advanced");
			}
			
			if(child->text(2)=="float")
			{
				param_->setValue(path,child->text(1).toDouble(),description,tags);
				String restrictions = child->data(2, Qt::UserRole).toString();
				vector<String> parts;
				if (restrictions.split(' ',parts))
				{
					if (parts[0]!="")
					{
						param_->setMinFloat(path,parts[0].toDouble());
					}
					if (parts[1]!="")
					{
						param_->setMaxFloat(path,parts[1].toDouble());
					}
				}
			}
			else if(child->text(2)=="string")
			{
				param_->setValue(path, child->text(1).toStdString(),description,tags);
				String restrictions = child->data(2, Qt::UserRole).toString();
				if(restrictions!="")
				{
					std::vector<String> parts;
					restrictions.split(',', parts);
					param_->setValidStrings(path,parts);
				}
			}
			else if(child->text(2) =="input file")
			{
				tags.push_back("input file");
				param_->setValue(path, child->text(1).toStdString(),description,tags);
				String restrictions = child->data(2, Qt::UserRole).toString();
				if(restrictions!="")
				{
					std::vector<String> parts;
					restrictions.split(',', parts);
					param_->setValidStrings(path,parts);
				}
			}
						else if(child->text(2) =="output file")
			{
				tags.push_back("output file");
				param_->setValue(path, child->text(1).toStdString(),description,tags);
				String restrictions = child->data(2, Qt::UserRole).toString();
				if(restrictions!="")
				{
					std::vector<String> parts;
					restrictions.split(',', parts);
					param_->setValidStrings(path,parts);
				}
			}
			else if(child->text(2)=="int")
			{
				param_->setValue(path, child->text(1).toInt(),description,tags);
				String restrictions = child->data(2, Qt::UserRole).toString();
				vector<String> parts;
				if (restrictions.split(' ',parts))
				{
					if (parts[0]!="")
					{
						param_->setMinInt(path,parts[0].toInt());
					}
					if (parts[1]!="")
					{
						param_->setMaxInt(path,parts[1].toInt());
					}
				}
			}
			String list;
			list = child->text(1).mid(1,child->text(1).length()-2); 
			StringList rlist = StringList::create(list);
			for(UInt i = 0; i < rlist.size(); ++i)
			{
				rlist[i] = rlist[i].trim();
			}
			if(child->text(2)=="string list")
			{
				param_->setValue(path,rlist,description,tags);
				String restrictions = child->data(2,Qt::UserRole).toString();
				if(restrictions!="")
				{
					vector<String> parts;
					restrictions.split(',',parts);
					param_->setValidStrings(path,parts);
				}
			}
			else if(child->text(2)=="input file list")
			{
				tags.push_back("input file");
				param_->setValue(path,rlist,description,tags);
				String restrictions = child->data(2,Qt::UserRole).toString();
				if(restrictions!="")
				{
					vector<String> parts;
					restrictions.split(',',parts);
					param_->setValidStrings(path,parts);
				}
			}
			else if(child->text(2)=="output file list")
			{
				tags.push_back("output file");
				param_->setValue(path,rlist,description,tags);
				String restrictions = child->data(2,Qt::UserRole).toString();
				if(restrictions!="")
				{
					vector<String> parts;
					restrictions.split(',',parts);
					param_->setValidStrings(path,parts);
				}
			}
			else if(child->text(2) =="double list")
			{
				param_->setValue(path,DoubleList::create(rlist),description,tags);
				String restrictions = child->data(2,Qt::UserRole).toString();
				vector<String> parts;
				if(restrictions.split(' ',parts))
				{
					if(parts[0]!= "")
					{
						param_->setMinFloat(path,parts[0].toFloat());
					}
					if(parts[1] != "")
					{
						param_->setMaxFloat(path,parts[1].toFloat());
					}
				}
			}
			else if(child->text(2) == "int list")
			{
				param_->setValue(path,IntList::create(rlist),description,tags);
				String restrictions = child->data(2,Qt::UserRole).toString();
				vector<String> parts;
				if(restrictions.split(' ',parts))
				{
					if(parts[0]!= "")
					{
						param_->setMinInt(path,parts[0].toInt());
					}
					if(parts[1] != "")
					{
						param_->setMaxInt(path,parts[1].toInt());
					}
				}			
			}
			
			// set description node description if the prefix matches
			for (map<String,String>::const_iterator it = section_descriptions.begin(); it!=section_descriptions.end(); ++it)
			{
				if (path.hasPrefix(it->first))
				{
					param_->setSectionDescription(it->first, it->second);
				}
			}
			section_descriptions.clear();
		}
		
		for (Int i = 0; i < child->childCount();++i)
		{
			storeRecursive_(child->child(i),path,section_descriptions);	//whole tree recursively
		}	
	}
	
	void ParamEditor::setModified(bool is_modified)
	{
		if (is_modified != modified_)
		{
			modified_ = is_modified;
			emit modified(modified_);
		}
	}

	bool ParamEditor::isModified() const
	{
		return modified_;
	}

	void ParamEditor::toggleAdvancedMode(bool advanced)
	{
		advanced_mode_ = advanced;
		
		stack<QTreeWidgetItem*> stack, node_stack;
		
		//show/hide items
		stack.push(tree_->invisibleRootItem());
		while(!stack.empty())
		{
			QTreeWidgetItem* current = stack.top();
			stack.pop();
			
			Int type = current->data(0,Qt::UserRole).toInt();
			if (type!=NODE) //ITEM
			{
				if (advanced_mode_ && type==ADVANCED_ITEM) //advanced mode
				{
					current->setHidden(false);
				}
				else if (!advanced_mode_ && type==ADVANCED_ITEM) //Normal mode
				{
					current->setHidden(true);
				}
			}
			else //NODE
			{
				for (Int i=0; i<current->childCount(); ++i)
				{
					stack.push(current->child(i));
				}
				
				if (advanced_mode_)
				{
					current->setHidden(false); //show all nodes in advanced mode
				}
				else
				{
					node_stack.push(current); //store node pointers in normal mode
				}
			}
		}
		
		//hide sections that have no visible items in normal mode
		while(!node_stack.empty())
		{
			QTreeWidgetItem* current = node_stack.top();
			node_stack.pop();
			
			bool has_visible_children = false;
			for (Int i=0; i<current->childCount(); ++i)
			{
				if (!current->child(i)->isHidden())
				{
					has_visible_children = true;
					break;
				}
			}
			if (!has_visible_children)
			{
				current->setHidden(true);
			}
		}
		
		//resize columns
		tree_->resizeColumnToContents(0);
		tree_->resizeColumnToContents(1);
		tree_->resizeColumnToContents(2);
		tree_->resizeColumnToContents(3);
	}

} // namespace OpenMS
