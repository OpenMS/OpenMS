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

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/SYSTEM/File.h>
#include <QtCore/QStringList>
#include <QtGui/QPushButton>
#include <QtGui/QComboBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QGridLayout>
#include <QtGui/QLabel>
#include <QtGui/QMessageBox>
#include <QtGui/QRadioButton>
#include <QtGui/QFileDialog>
#include <QtGui/QCheckBox>

using namespace std;

namespace OpenMS
{

	ToolsDialog::ToolsDialog( QWidget * parent, String tmp_dir, String default_dir, const LayerData* layer)
		: QDialog(parent),
			tmp_dir_(tmp_dir),
			default_dir_(default_dir)			
	{
		QGridLayout *main_grid=new QGridLayout(this);
		QLabel *label=NULL;

		label=new QLabel("TOPP tool:");
		main_grid->addWidget(label,0,0);
		QStringList list;
				
		if (layer->type==LayerData::DT_PEAK)
		{
			list<<"FileFilter"<<"FileInfo"<<"NoiseFilter"<<"BaselineFilter"<<"PeakPicker"<<"Resampler"<<"SpectraFilter"<<"MapNormalizer"<<"InternalCalibration"<<"ExternalCalibration"<<"FeatureFinder";
		}
		else if (layer->type==LayerData::DT_FEATURE)
		{
			list<<"FileFilter"<<"FileConverter"<<"FileInfo"<<"LabeledMatcher"<<"Decharger";
		}
		else if (layer->type==LayerData::DT_FEATURE_PAIR)
		{
			list<<"FileConverter"<<"LabeledMatcher"<<"Decharger"<<"FeaturePairSplitter";
		}
		list.push_front("<select>");
		tools_combo_=new QComboBox;
		tools_combo_->addItems(list);
		connect(tools_combo_,SIGNAL(activated(int)),this,SLOT(setTool_(int)));
		main_grid->addWidget(tools_combo_,0,1);
		
		label=new QLabel("input argument:");
		main_grid->addWidget(label,1,0);
		input_combo_=new QComboBox;
		input_combo_->setEnabled(false);
		main_grid->addWidget(input_combo_,1,1);
		
		label=new QLabel("output argument:");
		main_grid->addWidget(label,2,0);
		output_combo_=new QComboBox;
		output_combo_->setEnabled(false);
		main_grid->addWidget(output_combo_,2,1);
		
		
		QGridLayout* radio_grid = new QGridLayout;
		label=new QLabel("Output action:");
		output_radio_=new QRadioButton("None");
		radio_grid->addWidget(output_radio_,0,0);
		window_radio_=new QRadioButton("Open as new Window");
		radio_grid->addWidget(window_radio_,1,0);
		layer_radio_=new QRadioButton("Open as new Layer");
		layer_radio_->setChecked(true);
		radio_grid->addWidget(layer_radio_,2,0);
		main_grid->addWidget(label,3,0);
		main_grid->addLayout(radio_grid,3,1);

		//Add advanced mode check box		
		editor_=new ParamEditor(this);
		editor_->createShortcuts();
		main_grid->addWidget(editor_,4,0,1,4);
		QCheckBox* advanced = new QCheckBox("Show advanced parameters",this);
		main_grid->addWidget(advanced,5,3);
		connect(advanced,SIGNAL(toggled(bool)),editor_,SLOT(toggleAdvancedMode(bool)));
		main_grid->setColumnStretch(2,2);
		
		
		QHBoxLayout* hbox = new QHBoxLayout;
		QPushButton* load_button=new QPushButton(tr("&Load"));
		connect(load_button,SIGNAL(clicked()),this,SLOT(loadIni()));
		hbox->addWidget(load_button);
		QPushButton* store_button=new QPushButton(tr("&Store"));
		connect(store_button,SIGNAL(clicked()),this,SLOT(storeIni()));
		hbox->addWidget(store_button);
		hbox->addStretch();
		
		ok_button_= new QPushButton(tr("&Ok"));
		connect(ok_button_, SIGNAL(clicked()),this,SLOT(ok_()));
		hbox->addWidget(ok_button_);
		ok_button_->setEnabled(false);
		
		QPushButton* cancel_button=new QPushButton(tr("&Cancel"));
		connect(cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
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
		
		
		if(i==0)
		{
			
			ok_button_->setEnabled(false);
			input_combo_->setCurrentIndex(0);
			input_combo_->setEnabled(false);
			output_combo_->setCurrentIndex(0);
			output_combo_->setEnabled(false);
			output_combo_->setEnabled(false);
			editor_->clear();
			return;
		}
		
		String call = ToolsDialog::getTool()+" -write_ini "+tmp_dir_+"/in.ini -log "+tmp_dir_+"/ToolsDialog.log";
		
		if(system(call.c_str())!=0)
		{
			QMessageBox::critical(this,"Error",(String("Could not execute '")+call+"'!\nMake sure the TOPP tools are in your $PATH variable!").c_str());
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
				vis_param_.clear();
				editor_->clear();
				arg_map_.clear();
			}
			
			arg_param_.load((tmp_dir_+"/in.ini").c_str());
			
			vis_param_=arg_param_.copy(getTool()+":1:",true);
			
			editor_->loadEditable(vis_param_);
			
			String str;
			QStringList arg_list;
			for (Param::ParamIterator iter=arg_param_.begin();iter!=arg_param_.end();++iter)
			{
				str=iter.getName().substr(iter.getName().rfind("1:")+2,iter.getName().size());
				if(str.size()!=0 && str.find(":")==String::npos)
				{
					arg_map_.insert(make_pair(str,iter.getName()));
					arg_list<<QStringList(str.c_str());
				}
			}
			
			arg_list.push_front("<select>");
			input_combo_->clear();
			output_combo_->clear();
			input_combo_->addItems(arg_list);
			Int pos = arg_list.indexOf("in");
			if (pos!=-1)
			{
				input_combo_->setCurrentIndex(pos);
			}
			output_combo_->addItems(arg_list);
			pos = arg_list.indexOf("out");
			if (pos!=-1)
			{
				output_combo_->setCurrentIndex(pos);
			}
			editor_->setFocus(Qt::MouseFocusReason);
		}
	}
	
	void ToolsDialog::ok_()
	{
		if ((input_combo_->currentText()=="<select>" || output_combo_->currentText()=="<select>" || tools_combo_->currentText()=="<select>") && !noOutputAction())
		{
			QMessageBox::critical(this,"Error","You have to select a tool, an input argument and an output argument!");
		}
		else if((input_combo_->currentText()=="<select>" || tools_combo_->currentText()=="<select>") && noOutputAction())
		{
			QMessageBox::critical(this,"Error","You have to select a tool and an input argument!");
		}
		else
		{
			editor_->store();
			arg_param_.insert(getTool()+":1:",vis_param_);
			if(!File::writable(tmp_dir_+"/in.ini"))
			{
				QMessageBox::critical(this,"Error",(String("Could not write to '")+tmp_dir_+"/in.ini'!").c_str());
			}
			arg_param_.store(tmp_dir_+"/in.ini");
			input_string_=input_combo_->currentText().toStdString();
			output_string_=output_combo_->currentText().toStdString();
	
			accept();
		}
	}
	
	bool ToolsDialog::loadIni()
	{
		QString string;
		filename_=QFileDialog::getOpenFileName(this,tr("Open ini file"),default_dir_.c_str(),tr("ini files (*.ini);; all files (*.*)"));
		
		if(!filename_.isEmpty())
		{
			ok_button_->setEnabled(true);
			input_combo_->setEnabled(true);
			output_combo_->setEnabled(true);
			if(!arg_param_.empty())
			{
				arg_param_.clear();
				vis_param_.clear();
				editor_->clear();
				arg_map_.clear();
			}
			
			arg_param_.load(filename_.toStdString());
			Param::ParamIterator iter=arg_param_.begin();
			String str;
			string=iter.getName().substr(0,iter.getName().find(":")).c_str();
			Int pos = tools_combo_->findText(string);
			if (pos!=-1)
			{
				tools_combo_->setCurrentIndex(pos);
			}
			vis_param_=arg_param_.copy(getTool()+":1:",true);
			
			editor_->loadEditable(vis_param_);
			
			QStringList arg_list;
			for (Param::ParamIterator iter=arg_param_.begin();iter!=arg_param_.end();++iter)
			{
				str=iter.getName().substr(iter.getName().rfind("1:")+2,iter.getName().size());
				if(str.size()!=0 && str.find(":")==String::npos)
				{
					arg_map_.insert(make_pair(str,iter.getName()));
					arg_list<<QStringList(str.c_str());
				}
			}
			
			arg_list.push_front("<select>");
			input_combo_->clear();
			output_combo_->clear();
			input_combo_->addItems(arg_list);
			pos = arg_list.indexOf("in");
			if (pos!=-1)
			{
				input_combo_->setCurrentIndex(pos);
			}
			output_combo_->addItems(arg_list);
			pos = arg_list.indexOf("out");
			if (pos!=-1)
			{
				output_combo_->setCurrentIndex(pos);
			}
			
			string=QString("%1 - TOPP tools").arg(filename_);
			setWindowTitle(string.remove(0,string.lastIndexOf('/')+1));
			return true;
		}
		else return false;
	}
	
	bool ToolsDialog::storeIni()
	{
		filename_=QFileDialog::getSaveFileName(this,tr("Save ini file"),default_dir_.c_str(),tr("ini files (*.ini)"));
		if(!filename_.isEmpty() && !arg_param_.empty())
		{
			if(!filename_.endsWith(".ini")) filename_.append(".ini");
			editor_->store();
			arg_param_.insert(getTool()+":1:",vis_param_);
			if(!File::writable(filename_.toStdString()))
			{
				QMessageBox::critical(this,"Error",(String("Could not write to '")+filename_.toStdString()).c_str());
			}
			arg_param_.store(filename_.toStdString());
			QString str=QString("%1 - TOPP tools").arg(filename_);
			setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
			return true;
		}
		return false;
	}
	
	String ToolsDialog::getOutput()
	{
		return output_string_;
	}
	
	String ToolsDialog::getInput()
	{
		return input_string_;
	}
	
	bool ToolsDialog::openAsWindow()
	{
		return window_radio_->isChecked();
	}
	
	bool ToolsDialog::openAsLayer()
	{
		return layer_radio_->isChecked();
	}
	
	bool ToolsDialog::noOutputAction()
	{
		return output_radio_->isChecked();
	}
	
	String ToolsDialog::getTool()
	{
		return tools_combo_->currentText().toStdString();
	}
	
}
	
