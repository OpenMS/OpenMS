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
// $Maintainer: Stefan Rink $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/SYSTEM/File.h>
#include <QtCore/QStringList>
#include <QtGui/QPushButton>
#include <QtGui/QComboBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QGridLayout>
#include <QtGui/QLabel>
#include <QtGui/QMessageBox>
#include <QtGui/QRadioButton>


using namespace std;

namespace OpenMS
{

	ToolsDialog::ToolsDialog( QWidget * parent, String tmp_dir)
		: QDialog(parent),
			tmp_dir_(tmp_dir)
	{
		QHBoxLayout *hbox=NULL;
		QGridLayout *main_grid=NULL;
		QGridLayout *radio_grid=NULL;
		QLabel *label=NULL;
		
		QStringList list = TOPPBase::registerTools();
		list.push_front("<select>");
		tools_combo_=new QComboBox;
		tools_combo_->addItems(list);
		connect(tools_combo_,SIGNAL(activated(int)),this,SLOT(setTool_(int)));
		main_grid=new QGridLayout;
		main_grid->addWidget(tools_combo_,0,0,1,2);
		
		label=new QLabel("input file:");
		label->setAlignment(Qt::AlignHCenter|Qt::AlignVCenter);
		input_combo_=new QComboBox;
		input_combo_->setEnabled(false);
		main_grid->addWidget(label,1,0);
		main_grid->addWidget(input_combo_,1,1,1,1);
		
		label=new QLabel("output file:");
		label->setAlignment(Qt::AlignHCenter|Qt::AlignVCenter);
		output_combo_=new QComboBox;
		output_combo_->setEnabled(false);
		main_grid->addWidget(label,2,0);
		main_grid->addWidget(output_combo_,2,1,1,1);
		
		
		radio_grid=new QGridLayout;
		label=new QLabel("Open As:");
		label->setAlignment(Qt::AlignHCenter|Qt::AlignVCenter);
		window_radio_=new QRadioButton("New Window");
		radio_grid->addWidget(window_radio_,0,0);
		layer_radio_=new QRadioButton("New Layer");
		layer_radio_->setChecked(true);
		radio_grid->addWidget(layer_radio_,1,0);
		main_grid->addWidget(label,3,0);
		main_grid->addLayout(radio_grid,3,1);
		
		
		
		
		editor_=new ParamEditor;
		
		main_grid->addWidget(editor_,4,0,2,4);
		
		
		hbox=new QHBoxLayout;
		hbox->addStretch();
		
		ok_button_= new QPushButton(tr("&Ok"));
		connect(ok_button_, SIGNAL(clicked()),this,SLOT(ok_()));
		hbox->addWidget(ok_button_);
		ok_button_->setEnabled(false);
		
		QPushButton* cancel_button=new QPushButton(tr("&Cancel"));
		connect(cancel_button,SIGNAL(clicked()),this,SLOT(cancel_()));
		hbox->addWidget(cancel_button);
		main_grid->addLayout(hbox,6,0,1,4);
		
		setLayout(main_grid);
		
		setWindowTitle(tr("TOPP tools"));
		
	}
	
	ToolsDialog::~ToolsDialog()
	{
	
	}
	
	
	
	void ToolsDialog::setTool_(int i)
	{
		if(i)
		{
			String call = ToolsDialog::getTool()+" -write_ini "+tmp_dir_+"/in.ini";
			if(system(call.c_str())!=0)
			{
				QMessageBox::critical(this,"Error",(String("Could not execute '")+call+"'!").c_str());
			}
			else if(!File::exists(tmp_dir_+"/in.ini"))
			{
				QMessageBox::critical(this,"Error",(String("Could not open '")+tmp_dir_+"/in.ini'!").c_str());
			}
			else
			{
				ok_button_->setEnabled(true);
				input_combo_->setEnabled(true);
				output_combo_->setEnabled(true);
				if(!arg_param_.empty())
				{
					arg_param_.clear();
					editor_->deleteAll();
					arg_map_.clear();
				}
				arg_param_.load((tmp_dir_+"/in.ini").c_str());
				editor_->loadEditable(arg_param_);
				String str;
				QStringList arg_list;
				for (Param::ConstIterator iter=arg_param_.begin();iter!=arg_param_.end();++iter)
				{
					str=iter->first.substr(iter->first.rfind("1:")+2,iter->first.size());
					if(str.size()!=0 && str.find(":")==String::npos)
					{
						arg_map_.insert(make_pair(str,iter->first));
						arg_list<<QStringList(str.c_str());
					}
				}
				arg_list.push_front("<select>");
				input_combo_->clear();
				output_combo_->clear();
				input_combo_->addItems(arg_list);
				output_combo_->addItems(arg_list);
			}

		}
	}
	
	void ToolsDialog::ok_()
	{
		editor_->store();
		arg_param_.store(tmp_dir_+"/in.ini");
		setInput_();
		setOutput_();
		
		accept();
	}
	
	void ToolsDialog::cancel_()
	{
		reject();
	}
	
	void ToolsDialog::setInput_()
	{
		if(!input_combo_->currentIndex() && input_combo_->currentIndex()!=-1)
		{
			input_string_=input_combo_->currentText().toStdString();
		}
		else if(!arg_param_.getValue(arg_map_["in"]).toString().empty())
		{
			input_string_=arg_param_.getValue(arg_map_["in"]).toString();
		}
		else
		{
			input_string_="in";
		}
	}
	
	void ToolsDialog::setOutput_()
	{
		if(!output_combo_->currentIndex() && output_combo_->currentIndex()!=-1)
		{
			output_string_=output_combo_->currentText().toStdString();
		}
		else if(!arg_param_.getValue(arg_map_["out"]).toString().empty())
		{
			output_string_=arg_param_.getValue(arg_map_["out"]).toString();
		}
		else
		{
			output_string_="out";
		}
	}
	
	String ToolsDialog::getOutput()
	{
		return output_string_;
	}
	
	String ToolsDialog::getInput()
	{
		return input_string_;
	}
	
	bool ToolsDialog::isWindow()
	{
		return window_radio_->isChecked();
	}
	
	String ToolsDialog::getTool()
	{
		return tools_combo_->currentText().toStdString();
	}
	
}
	
