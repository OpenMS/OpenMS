// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
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
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

using namespace std;

namespace OpenMS
{

	ToolsDialog::ToolsDialog( QWidget * parent, String ini_file, String default_dir, LayerData::DataType type)
		: QDialog(parent),
			ini_file_(ini_file),
			default_dir_(default_dir)			
	{
		QGridLayout *main_grid=new QGridLayout(this);
		QLabel *label=NULL;

		label=new QLabel("TOPP tool:");
		main_grid->addWidget(label,0,0);
		QStringList list;
				
		if (type==LayerData::DT_PEAK)
		{
			list<<"FileFilter"<<"FileInfo"<<"NoiseFilter"<<"BaselineFilter"<<"PeakPicker"<<"Resampler"<<"SpectraFilter"<<"MapNormalizer"<<"InternalCalibration"<<"TOFCalibration"<<"FeatureFinder"<<"SILACAnalyzer"<<"ITRAQAnalyzer";
		}
		else if (type==LayerData::DT_FEATURE)
		{
			list<<"FileFilter"<<"FileConverter"<<"FileInfo"<<"Decharger"<<"FeatureLinker";
		}
		else if (type==LayerData::DT_CONSENSUS)
		{
			list<<"FileFilter"<<"FileConverter"<<"FileInfo";
		}
		else if (type==LayerData::DT_CHROMATOGRAM)
		{
		  //TODO CHROM
		}
		//sort list alphabetically
		list.sort();
		list.push_front("<select tool>");
		tools_combo_=new QComboBox;
		tools_combo_->setMinimumWidth(150);
		tools_combo_->addItems(list);
		connect(tools_combo_,SIGNAL(activated(int)),this,SLOT(setTool_(int)));
		connect(tools_combo_,SIGNAL(currentIndexChanged(int)),this,SLOT(updateTypes_()));

		main_grid->addWidget(tools_combo_,0,1);

		//type combobox
		type_combo_=new QComboBox;
		type_combo_->setMinimumWidth(150);
		connect(type_combo_,SIGNAL(activated(int)),this,SLOT(setType_(int)));
		list.clear();
		list << "<select type>";
		type_combo_->addItems(list);

		type_combo_->setEnabled(false);
		main_grid->addWidget(type_combo_,0,2);
		
		
		label=new QLabel("input argument:");
		main_grid->addWidget(label,1,0);
		input_combo_=new QComboBox;
		main_grid->addWidget(input_combo_,1,1);
		
		label=new QLabel("output argument:");
		main_grid->addWidget(label,2,0);
		output_combo_=new QComboBox;
		main_grid->addWidget(output_combo_,2,1);

		//Add advanced mode check box		
		editor_=new ParamEditor(this);
		main_grid->addWidget(editor_,3,0,1,5);		
		
		QHBoxLayout* hbox = new QHBoxLayout;
		QPushButton* load_button=new QPushButton(tr("&Load"));
		connect(load_button,SIGNAL(clicked()),this,SLOT(loadINI_()));
		hbox->addWidget(load_button);
		QPushButton* store_button=new QPushButton(tr("&Store"));
		connect(store_button,SIGNAL(clicked()),this,SLOT(storeINI_()));
		hbox->addWidget(store_button);
		hbox->addStretch();
		
		ok_button_= new QPushButton(tr("&Ok"));
		connect(ok_button_, SIGNAL(clicked()),this,SLOT(ok_()));
		hbox->addWidget(ok_button_);
		
		QPushButton* cancel_button=new QPushButton(tr("&Cancel"));
		connect(cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
		hbox->addWidget(cancel_button);
		main_grid->addLayout(hbox,5,0,1,5);
		
		setLayout(main_grid);
		
		setWindowTitle(tr("TOPP tools"));
		disable_();		
	}
	
	ToolsDialog::~ToolsDialog()
	{
	
	}
	
	void ToolsDialog::setType_(int i)
	{	
		if (i==0)		//no type selected
		{
			disable_();
			editor_->clear();			
			return;
		}
		
		String call = getTool()+" -write_ini "+ini_file_+" -log "+ini_file_+".log";

		if (i!=-1)
		{
			call += " -type " + String(type_combo_->currentText());
		}
		
		if(system(call.c_str())!=0)
		{
			QMessageBox::critical(this,"Error",(String("Could not execute '")+call+"'!\n\nMake sure the TOPP tools are in your $PATH variable, that you have write permission in the temporary file path, and that there is space left in the temporary file path.").c_str());
		}
		else if(!File::exists(ini_file_))
		{
			QMessageBox::critical(this,"Error",(String("Could not open '")+ini_file_+"'!").c_str());
		}
		else
		{
			enable_();
			if(!arg_param_.empty())
			{
				arg_param_.clear();
				vis_param_.clear();
				editor_->clear();
				arg_map_.clear();
			}
			
			arg_param_.load((ini_file_).c_str());
			
			vis_param_=arg_param_.copy(getTool()+":1:",true);
			vis_param_.remove("log");
			vis_param_.remove("no_progress");
			vis_param_.remove("debug");
			
			//set selected type
			if (i!=-1)
			{
				vis_param_.setValue("type",String(type_combo_->currentText()));
			}
			
			editor_->load(vis_param_);
			
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
			if (pos!=-1 && getTool()!="FileInfo")
			{
				output_combo_->setCurrentIndex(pos);
			}
			editor_->setFocus(Qt::MouseFocusReason);
		}
	}
	
	void ToolsDialog::setTool_(int i)
	{
		editor_->clear();

		//no tool selected		
		if(i==0)
		{
			disable_();
			return;
		}

		if (ToolHandler::getTypes(getTool()).empty())
		{
			setType_(-1);
		}
	}

	void ToolsDialog::disable_()
	{
		ok_button_->setEnabled(false);
		input_combo_->setCurrentIndex(0);
		input_combo_->setEnabled(false);
		output_combo_->setCurrentIndex(0);
		output_combo_->setEnabled(false);
	}

	void ToolsDialog::enable_()
	{
		ok_button_->setEnabled(true);
		input_combo_->setEnabled(true);
		output_combo_->setEnabled(true);
	}

	void ToolsDialog::updateTypes_()
	{
		//clean up
		type_combo_->clear();
		
		//get type list
		StringList type_list = ToolHandler::getTypes(tools_combo_->currentText());
		
		//no types => disable combo box
		if (type_list.empty())
		{
			type_combo_->setEnabled(false);
			enable_();
		}
		//types => fill combo box
		else
		{
			QStringList list;
			list << "<select type>";		
			for (Size i=0; i<type_list.size(); ++i)
			{
				list << type_list[i].toQString();
			}

			type_combo_->addItems(list);
			type_combo_->setEnabled(true);
			disable_();
		}
	}


	void ToolsDialog::ok_()
	{
		if (input_combo_->currentText()=="<select>" || tools_combo_->currentText()=="<select>")
		{
			QMessageBox::critical(this,"Error","You have to select a tool and an input argument!");
		}
		else
		{
			editor_->store();
			arg_param_.insert(getTool()+":1:",vis_param_);
			if(!File::writable(ini_file_))
			{
				QMessageBox::critical(this,"Error",(String("Could not write to '")+ini_file_+"'!").c_str());
			}
			arg_param_.store(ini_file_);
			accept();
		}
	}
	
	void ToolsDialog::loadINI_()
	{
		QString string;
		filename_=QFileDialog::getOpenFileName(this,tr("Open ini file"),default_dir_.c_str(),tr("ini files (*.ini);; all files (*.*)"));
		//not file selected
		if(filename_.isEmpty())
		{
			return;
		}
		enable_();
		if(!arg_param_.empty())
		{
			arg_param_.clear();
			vis_param_.clear();
			editor_->clear();
			arg_map_.clear();
		}
		try
		{
			arg_param_.load(filename_.toStdString());
		}
		catch(Exception::BaseException e)
		{
			QMessageBox::critical(this,"Error",(String("Error loading INI file: ")+e.getMessage()).c_str());
			arg_param_.clear();
			return;
		}
		//set tool combo
		Param::ParamIterator iter=arg_param_.begin();
		String str;
		string=iter.getName().substr(0,iter.getName().find(":")).c_str();
		Int pos = tools_combo_->findText(string);
		if (pos==-1)
		{
			QMessageBox::critical(this,"Error",(String("Cannot apply '")+string+"' tool to this layer type. Aborting!").c_str());
			arg_param_.clear();
			return;
		}
		tools_combo_->setCurrentIndex(pos);
		//Extract the required parameters
		vis_param_=arg_param_.copy(getTool()+":1:",true);
		vis_param_.remove("log");
		vis_param_.remove("no_progress");
		vis_param_.remove("debug");
		//load data into editor
		editor_->load(vis_param_);			
		//special treatment for tools with type
		if (!ToolHandler::getTypes(string).empty())
		{
			String type;
			if (vis_param_.exists("type")) type = vis_param_.getValue("type");
			Int pos = type_combo_->findText(type.toQString());
			if (pos!=-1)
			{
				type_combo_->setCurrentIndex(pos);
				enable_();
			}
			else
			{
				type_combo_->setCurrentIndex(0);
				disable_();
				editor_->clear();
			}
		}
		
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
		if (pos!=-1 && getTool()!="FileInfo")
		{
			output_combo_->setCurrentIndex(pos);
		}
	}
	
	void ToolsDialog::storeINI_()
	{
		//nothing to save
		if (arg_param_.empty()) return;
		filename_=QFileDialog::getSaveFileName(this,tr("Save ini file"),default_dir_.c_str(),tr("ini files (*.ini)"));
		//not file selected
		if(filename_.isEmpty())
		{
			return;
		}
		if(!filename_.endsWith(".ini")) filename_.append(".ini");
		editor_->store();
		arg_param_.insert(getTool()+":1:",vis_param_);
		try
		{
			arg_param_.store(filename_.toStdString());
		}
		catch(Exception::BaseException e)
		{
			QMessageBox::critical(this,"Error",(String("Error storing INI file: ")+e.getMessage()).c_str());
			return;
		}
	}
	
	String ToolsDialog::getOutput()
	{
		if (output_combo_->currentText()=="<select>") return "";
		return output_combo_->currentText();
	}
	
	String ToolsDialog::getInput()
	{
		return input_combo_->currentText();
	}
	
	String ToolsDialog::getTool()
	{
		return tools_combo_->currentText().toStdString();
	}
	
}
	
