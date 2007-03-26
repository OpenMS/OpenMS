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
#include <QtCore/QString>
#include <QtCore/QStringList>
#include <QtGui/QPushButton>
#include <QtGui/QComboBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QVBoxLayout>
#include <QtGui/QLabel>
#include <QtGui/QMessageBox>
#include <QtGui/QRadioButton>
#include <QtGui/QLineEdit>
#include <QtCore/QDir>


using namespace std;

namespace OpenMS
{

	ToolsDialog::ToolsDialog( QWidget * parent, String tmp_dir)
		: QDialog(parent),
			tmp_dir_(tmp_dir)
	{
		QVBoxLayout *main=new QVBoxLayout;
		QVBoxLayout *vbox=NULL;
		QHBoxLayout *hbox=NULL;
		QGridLayout *grid=NULL;
		QValidator *validator=NULL;
		QLabel *label=NULL;
		
		QStringList list = TOPPBase::registerTools();
		list.push_front("<select>");
		tools_combo_=new QComboBox;
		tools_combo_->addItems(list);
		connect(tools_combo_,SIGNAL(activated(int)),this,SLOT(setTool_(int)));
		hbox=new QHBoxLayout;
		hbox->addWidget(tools_combo_);
		hbox->addStretch();
		main->addLayout(hbox);
		
		label=new QLabel("input file:");
		input_edit_=new QLineEdit;
		validator = new QRegExpValidator(QRegExp("\\w+"),this);
		input_edit_->setValidator(validator);
		grid=new QGridLayout;
		grid->addWidget(label,0,0);
		grid->addWidget(input_edit_,0,1);
		
		label=new QLabel("output file:");
		output_edit_=new QLineEdit;
		validator = new QRegExpValidator(QRegExp("\\w+"),this);
		output_edit_->setValidator(validator);
		grid->addWidget(label,1,0);
		grid->addWidget(output_edit_,1,1);
		main->addLayout(grid);
	
		label=new QLabel("Open As:");
		hbox=new QHBoxLayout;
		vbox=new QVBoxLayout;
		window_radio_=new QRadioButton("New Window");
		vbox->addWidget(window_radio_);
		layer_radio_=new QRadioButton("New Layer");
		layer_radio_->setChecked(true);
		vbox->addWidget(layer_radio_);
		hbox->addWidget(label);
		hbox->addLayout(vbox);
		hbox->addStretch();
		main->addLayout(hbox);
		
		
		
		
		editor_=new ParamEditor;
		
		main->addWidget(editor_);
		
		hbox=new QHBoxLayout;
		hbox->addStretch();
		
		ok_button_= new QPushButton(tr("&Ok"));
		connect(ok_button_, SIGNAL(clicked()),this,SLOT(ok_()));
		hbox->addWidget(ok_button_);
		ok_button_->setEnabled(false);
		
		QPushButton* cancel_button=new QPushButton(tr("&Cancel"));
		connect(cancel_button,SIGNAL(clicked()),this,SLOT(cancel_()));
		hbox->addWidget(cancel_button);
		
		main->addLayout(hbox);
		
		setLayout(main);
		
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
				if(!arg_param_.empty())
				{
					arg_param_.clear();
					editor_->deleteAll();
					arg_map_.clear();
				}
				arg_param_.load((tmp_dir_+"/in.ini").c_str());
				editor_->loadEditable(arg_param_);
				String str;
				for (Param::ConstIterator iter=arg_param_.begin();iter!=arg_param_.end();++iter)
				{
					str=iter->first.substr(iter->first.rfind("1:")+2,iter->first.size());
					if(str.size()!=0 && str.find(":")==String::npos)
					{
						arg_map_.insert(make_pair(str,iter->first));
					}
				}
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
		if(!input_edit_->text().isEmpty())
		{
			input_string_=input_edit_->text().toStdString();
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
		if(!output_edit_->text().isEmpty())
		{
			output_string_=output_edit_->text().toStdString();
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
	
