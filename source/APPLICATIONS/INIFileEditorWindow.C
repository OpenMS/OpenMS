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

#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/APPLICATIONS/INIFileEditorWindow.h>

#include <QtGui/QApplication>
#include <QtGui/QMainWindow>
#include <QtGui/QToolBar>
#include <QtCore/QString>
#include <QtGui/QFileDialog>
#include <QtGui/QIcon>

#include "../VISUAL/ICONS/save.xpm"
#include "../VISUAL/ICONS/saveas.xpm"
#include "../VISUAL/ICONS/open.xpm"
#include "../VISUAL/ICONS/TOPPView.xpm"

using namespace OpenMS;


		 INIFileEditorWindow::INIFileEditorWindow(QWidget *parent) : QMainWindow(parent)
		{
			
			
			setWindowTitle("INIFileEditor");
			setWindowIcon(QIcon(XPM_toppview));
			editor_=new ParamEditor;
			setCentralWidget(editor_);
			
			
			
			toolbar_=addToolBar(tr("&File"));
			toolbar_->addAction(QIcon(open_xpm),tr("&Open ini file"),this,SLOT(openFile_()));
			toolbar_->addAction(QIcon(save_xpm),tr("&Save ini file"),this,SLOT(saveFile_()));
			toolbar_->addAction(QIcon(saveas_xpm),tr("&Save as"),this,SLOT(saveFileAs_()));
		}

		bool INIFileEditorWindow::openFile_()
		{
			filename_=QFileDialog::getOpenFileName(this,tr("Open ini file"),".",tr("ini files (*.ini)"));
			if(!filename_.isEmpty())
			{
				param_.load(filename_.toStdString());
				editor_->loadEditable(param_);
				QString str=QString("%1 - INIFileEditor").arg(filename_);
				setWindowTitle(str.remove(0,str.lastIndexOf('/')+1));
				return true;
			}
			else return false;
		}
		bool INIFileEditorWindow::saveFile_()
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
				return saveFileAs_();
			}
		}
	
		bool INIFileEditorWindow::saveFileAs_()
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
