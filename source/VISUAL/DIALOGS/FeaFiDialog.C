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
#include <OpenMS/VISUAL/DIALOGS/FeaFiDialog.h>
#include <sstream>
#include <QtGui/QFileDialog>
#include <QtGui/QPushButton>
#include <QtGui/QLabel>

using namespace std;

namespace OpenMS
{

	FeaFiDialog::FeaFiDialog( QWidget * parent)
		: QDialog(parent), 
			finder_()
	{
		setupUi(this);
		
		resize(sizeHint());
		start_button->setEnabled(false);
	}
	
	FeaFiDialog::~FeaFiDialog()
	{
	
	}
	
	FeatureFinder& FeaFiDialog::getFeatureFinder()
	{
		return finder_;
	}
	
	
	void FeaFiDialog::loadParamFile()
	{
		QString dir = QDir::current().path();
		dir = dir.left(dir.indexOf("OpenMS")).append("OpenMS/source/APPLICATIONS/FEATUREFINDER/");
	 	QString file = QFileDialog::getOpenFileName(this, "Open file", dir, "Parameters (*.ini)");
	
		if (!file.isEmpty() && file!="")
		{
			Param param;
			param.load(file.toAscii().data());
			bool isParamValid = finder_.setParam(param);
			if (isParamValid)
			{
				std::ostringstream label;
				label << finder_;
				param_label->setText(label.str().c_str());
				param_label->setAlignment( Qt::AlignVCenter|Qt::AlignLeft );
				start_button->setEnabled(true);
				start_button->setFocus();
			}
			else
			{
				param_label->setText(tr("Chosen param file was not valid.\nNo Feature-Finder loaded yet."));
			}	
		}
	}

} //namespace OpenMS
