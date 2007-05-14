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
		: QMainWindow(parent)
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
		
		edit_ = new QMenu("&Edit",this);
		menuBar()->addMenu(edit_);
		edit_->addAction("&Copy",editor_,SLOT(copySubTree()), Qt::CTRL+Qt::Key_C);
		edit_->addAction("&Cut",editor_,SLOT(cutSubTree()), Qt::CTRL+Qt::Key_X);
		edit_->addAction("&Paste",editor_,SLOT(pasteSubTree()), Qt::CTRL+Qt::Key_V);
		edit_->addAction("&Delete",editor_,SLOT(deleteItem()), Qt::Key_Delete);
		edit_->addSeparator();
		edit_->addAction("Insert new &Section",editor_,SLOT(insertNode()), Qt::Key_S);
		edit_->addAction("Insert new &Value",editor_,SLOT(insertItem()), Qt::Key_V);
		//edit_->setEnabled(false);
		
		connect(editor_,SIGNAL(modified(bool)),this,SLOT(updateWindowTitle(bool)));
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
		if(filename_.isEmpty())
		{
			QMessageBox::warning(this,"No ini-file!","You have to open an ini-file before saving!");
			return false;
		}
		else if(editor_->isNameEmpty())
		{
			QMessageBox::warning(this,"Empty Name","You have to enter a name before saving!");
			return false;
		}
		
		editor_->store();
		param_.store(filename_.toStdString());
		QString str=QString("%1 - INIFileEditor").arg(filename_);
		setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
		
		return true;
	}
	
	
	bool INIFileEditorWindow::saveFileAs()
	{
		filename_=QFileDialog::getSaveFileName(this,tr("Save ini file"),".",tr("ini files (*.ini)"));
		if(!filename_.isEmpty() && !editor_->isNameEmpty())
		{
			if(!filename_.endsWith(".ini")) filename_.append(".ini");
			editor_->store();
			param_.store(filename_.toStdString());
			QString str=QString("%1 - INIFileEditor").arg(filename_);
			setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
			return true;
		}
		else if(editor_->isNameEmpty())
		{
			QMessageBox::warning(this,"Empty Name","You have to enter a name before saving!");
		}
		return false;
	}
	
	void INIFileEditorWindow::closeEvent(QCloseEvent* /*event*/)
	{
		if(editor_->isModified())
		{
			if (QMessageBox::question(this,"Save?","Do you want to save your changes?",QMessageBox::Ok|QMessageBox::Cancel)==QMessageBox::Ok)
			{
				saveFile();
			}
		}
	}
	void INIFileEditorWindow::updateWindowTitle(bool update)
	{
		if(update)
		{
			QString str=QString("%1 * - INIFileEditor").arg(filename_);
			setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
		}
		else
		{
			QString str=QString("%1 - INIFileEditor").arg(filename_);
			setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
		}
	}
}

