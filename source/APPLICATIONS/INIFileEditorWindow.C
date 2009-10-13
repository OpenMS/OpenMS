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

#include <OpenMS/APPLICATIONS/INIFileEditorWindow.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtGui/QToolBar>
#include <QtCore/QString>
#include <QtGui/QFileDialog>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QMessageBox>
#include <QtGui/QCloseEvent>
#include <QtGui/QGridLayout>
#include <QtGui/QCheckBox>

using namespace std;

namespace OpenMS
{

	INIFileEditorWindow::INIFileEditorWindow(QWidget *parent) 
		: QMainWindow(parent),
			current_path_(".")
	{
		setWindowTitle("INIFileEditor");
	  setWindowIcon(QIcon(":/INIFileEditor.png"));

		//create central widget and layout
		QWidget* central_widget = new QWidget;
		setCentralWidget(central_widget);
		QGridLayout* layout = new QGridLayout(central_widget);
		
		//create advanced check box and ParamEditor and connect them
		editor_=new ParamEditor(central_widget);
		layout->addWidget(editor_,0,0,1,2);
		
		QMenu* file = new QMenu("&File",this);
		menuBar()->addMenu(file);
		file->addAction("&Open",this,SLOT(openFile()), Qt::CTRL+Qt::Key_O);
		file->addSeparator();
		file->addAction("&Save",this,SLOT(saveFile()), Qt::CTRL+Qt::Key_S);
		file->addAction("Save &As",this,SLOT(saveFileAs()));
		file->addSeparator();
		file->addAction("&Quit",this,SLOT(close()));
		
		// we connect the "changes state"(changes made/no changes) signal from the ParamEditor to the window title updating slot
		connect(editor_,SIGNAL(modified(bool)),this,SLOT(updateWindowTitle(bool)));
		
		setMinimumSize(600,600);
	}
	
	bool INIFileEditorWindow::openFile(const String& filename)
	{
		if (filename=="")
		{
			filename_=QFileDialog::getOpenFileName(this,tr("Open ini file"),current_path_.toQString(),tr("ini files (*.ini);; all files (*.*)"));
		}
		else
		{
			filename_ = filename.c_str();
		}
		
		if(!filename_.isEmpty())
		{
			if (File::readable(filename_.toStdString()))
			{
				param_.clear();
				param_.load(filename_.toStdString());
				editor_->load(param_);
				updateWindowTitle(editor_->isModified());
				return true;
			}
			else
			{
				QMessageBox::critical(this,"Error opening file",("The file '"+filename_.toStdString()+"' does not exist or is not readable!").c_str());		
			}
		}
		return false;
	}
	
	
	bool INIFileEditorWindow::saveFile()
	{
		if(filename_.isEmpty())
		{
			return false;
		}
		
		editor_->store();

		param_.store(filename_.toStdString());
		updateWindowTitle(editor_->isModified());
		return true;
	}
	
	
	bool INIFileEditorWindow::saveFileAs()
	{
		filename_=QFileDialog::getSaveFileName(this,tr("Save ini file"),current_path_.toQString(),tr("ini files (*.ini)"));
		if(!filename_.isEmpty())
		{
			if(!filename_.endsWith(".ini")) filename_.append(".ini");
			
			editor_->store();
			
			param_.store(filename_.toStdString());
			updateWindowTitle(editor_->isModified());
			return true;
		}
		return false;
	}
	
	void INIFileEditorWindow::closeEvent(QCloseEvent* event)
	{
		if(editor_->isModified())
		{
			QMessageBox::StandardButton result=QMessageBox::question(this,"Save?","Do you want to save your changes?",QMessageBox::Ok|QMessageBox::Cancel|QMessageBox::Discard);
			if (result==QMessageBox::Ok)
			{
				if(saveFile())
				{
					event->accept();
				}
				else
				{
					event->ignore();
				}
			}
			else if(result==QMessageBox::Cancel)
			{
				event->ignore();
			}
			else
			{
				event->accept();
			}
		}
		else
		{
			event->accept();
		}
	}
	
	void INIFileEditorWindow::updateWindowTitle(bool update)
	{
		//update window title			
		if(update)
		{
			setWindowTitle((File::basename(filename_) + " * - INIFileEditor").toQString());
		}
		else
		{
			setWindowTitle((File::basename(filename_) + " - INIFileEditor").toQString());
		}
		
		//update last path as well
		current_path_ = File::path(filename_);
	}
}

