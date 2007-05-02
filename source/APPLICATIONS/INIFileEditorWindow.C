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

#include <OpenMS/APPLICATIONS/INIFileEditorWindow.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/SYSTEM/File.h>
#include <QtGui/QToolBar>
#include <QtCore/QString>
#include <QtGui/QFileDialog>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QMessageBox>
#include <QtGui/QCloseEvent>

namespace OpenMS
{

	INIFileEditorWindow::INIFileEditorWindow(QWidget *parent) 
		: QMainWindow(parent),changed_(false)
	{
		setWindowTitle("INIFileEditor");
		editor_=new ParamEditor;
		setCentralWidget(editor_);
		
		QMenu* file = new QMenu("&File",this);
		menuBar()->addMenu(file);
		file->addAction("&Open",this,SLOT(openFile()), Qt::CTRL+Qt::Key_O);
		file->addSeparator();
		file->addAction("&Save",this,SLOT(saveFile()), Qt::CTRL+Qt::Key_S);
		file->addAction("Save &As",this,SLOT(saveFileAs()));
		file->addSeparator();
		file->addAction("&Quit",this,SLOT(close()));
		connect(editor_,SIGNAL(itemChanged ( QTreeWidgetItem *, int)),this,SLOT(setChanged(QTreeWidgetItem *, int)));
	}
	
	
	bool INIFileEditorWindow::openFile(const String& filename)
	{
		if (filename=="")
		{
			filename_=QFileDialog::getOpenFileName(this,tr("Open ini file"),".",tr("ini files (*.ini);; all files (*.*)"));
		}
		else
		{
			filename_ = filename.c_str();
		}
		
		if(!filename_.isEmpty())
		{
			if (File::readable(filename_.toStdString()))
			{
				param_.load(filename_.toStdString());
				editor_->loadEditable(param_);
				QString str=QString("%1 - INIFileEditor").arg(filename_);
				setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
				return true;
			}
			else
			{
				QMessageBox::critical(this,"Error opeing file",("The file '"+filename_.toStdString()+"' does not exist or is not readable!").c_str());		
			}
		}
		return false;
	}
	
	
	bool INIFileEditorWindow::saveFile()
	{
		if(!filename_.isEmpty())
		{
			editor_->store();
			param_.store(filename_.toStdString());
			QString str=QString("%1 - INIFileEditor").arg(filename_);
			setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
			
			return true;
		}
		else
		{
			return saveFileAs();
		}
	}
	
	
	bool INIFileEditorWindow::saveFileAs()
	{
		filename_=QFileDialog::getSaveFileName(this,tr("Save ini file"),".",tr("ini files (*.ini)"));
		if(!filename_.isEmpty())
		{
			if(!filename_.endsWith(".ini")) filename_.append(".ini");
			editor_->store();
			param_.store(filename_.toStdString());
			QString str=QString("%1 - INIFileEditor").arg(filename_);
			setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
			return true;
		}
		return false;
	}
	
	void INIFileEditorWindow::closeEvent(QCloseEvent* /*event*/)
	{
		if(isChanged())
		{
			if (QMessageBox::question(this,"Save?","Do you want to save your changes?",QMessageBox::Ok|QMessageBox::Cancel)==QMessageBox::Ok)
			{
				saveFile();
			}
		}
	}
	void INIFileEditorWindow::setChanged(QTreeWidgetItem * item, int column)
	{
		if(item->data(column,33).isValid() && item->data(column,33)!=item->data(column,Qt::DisplayRole))
		{
			changed_=true;
			QString str=QString("%1 * - INIFileEditor").arg(filename_);
			setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
		}
		else
		{
			changed_=false;
			QString str=QString("%1 - INIFileEditor").arg(filename_);
			setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
		}
	}
	
	bool INIFileEditorWindow::isChanged()
	{
		return changed_;
	}
	
}

